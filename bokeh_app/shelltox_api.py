	
# coding: utf-8

# # Shellfish Toxicity Modeling api
# This Jupyter Notebook has a number of functions useful for plotting shellfish toxicity and generating a shellfish toxicity ~~forecast~~ model.

# 
# ### The Model
# 
#     T(t) = c1 * T(t-1) + c2 * C(t-τ)
# where
# * T(t) is the modeled toxicity for day t. AKA uptake. 
# * T(0) by default is the first available measured toxicity value (where there is ALSO available ESP data) is used. It can however be set to any desired starting value.
# * C(t) is an ESP cell count for day t, values are interpolated if missing (minus tau τ days to account for drift lag, if any) 
# * c1 (aka 1-γ) is the toxicity constant, gamma γ is Depuration constant
# * c2 (aka α) is the cell count constant, as α*C() this is our uptake
# 
# c1 and c2 are found by performing a gradient descent search over an option space of thousands possible values before finally selecting the pair resulting in the least root mean square deviation between the available measured toxicity values and modeled toxicity values for the same dates as the available measured values.

#  

# In[1]:


## IMPORTS ##
import pandas as pd
from math import isnan
from pprint import pprint
from statistics import mean
from math import sqrt
from io import StringIO

## BOKEH IMPORTS & SETUP ##
from bokeh.plotting import figure, output_notebook, show, curdoc
from bokeh.layouts import column, row, WidgetBox
from bokeh.models import ColumnDataSource, CustomJS
from bokeh.models import LinearAxis, Range1d, DatetimeTickFormatter
from bokeh.models import Span, Slope, DateSlider, WheelZoomTool, PreText
from bokeh.models import Panel, Tabs, DataTable, TableColumn, DateFormatter, NumberFormatter
from bokeh.palettes import Colorblind
COLORS = Colorblind[7]
COLORS.append(COLORS.pop(2))
COLORS.insert(1,COLORS.pop(2))

#output_notebook()  # necessary for bokeh plots to appear in notebook


# In[2]:


def load_dataset(fname):
    '''Loads a csv with date, ESP, TOX and optionally MODEL data columns'''
    
    # Ingest the csv, interpretting the first column as a dates index
    try: df = pd.read_csv(fname, index_col=0, parse_dates=True)
    except: df = pd.read_csv('data/{}'.format(fname), index_col=0, parse_dates=True)
    
    # for blank lines, the dates are interpretted as NaN which messes up later functions, so remove them
    df = df[df.index.notnull()]

    #removing duplicate dates
    if any(df.index.duplicated()):
        print('DF CONTAINED DUPLICATE ROWS')
        df = df[~df.index.duplicated(keep='first')]
    
    # insert rows where there are missing days to get a monotonic timeseries
    df = df.reindex( pd.date_range(start=df.index[0],end=df.index[-1]+pd.DateOffset(days=1)) )
    
    df.index.name = 'date'

    return df


# In[3]:


# TEST_DF 
# a simple artificial dataset to test certain functionalities

#df_test = pd.DataFrame({'model_vals':[pd.np.nan,pd.np.nan,150,240,123,134,200,100],
#                        'esp_vals':[pd.np.nan,222,333,444,444,555,333,444],
#                        'tox_vals':[220,33,44,33,22,44,33,220]})
#df_test.index = pd.DatetimeIndex(['2014-03-31', '2014-04-01', '2014-04-02', '2014-04-03', 
#                                  '2014-04-04', '2014-04-05', '2014-04-06', '2014-04-07'])
#df_test.index.name = 'date'


# In[4]:


# The primary plotting function

def plot_shelltox(df , title=None, width=800, height=450, show_plot=False, show_slider=True, x_start=None, x_end=None, hide_toolbar=False):
    '''
    Plots the columns in df against the date index of df
    "ESP" columns are formatted with solid lines and cross datapoints on the left-hand-side cell-count axis
    "MODEL" columns are formatted with dashed lines and dotted datapoints
    "og check TOX" and "og check MODEL" are formatted as big dot datapoints
    Any other columns (ie TOX shellfish sites) will get plotted with solid lines and dotted data points
    All but "ESP" columns are plotted against the toxicity axis
    
    The plot also features a horizontal dotted red line at 80 µg/100g toxicity. 
    Beyond this toxicity, fisheries close their shellfish beds. 
    
    title changes the title of the plot, the default is "Toxicity and Cell Count"
    width and height determines the dimensions of the plot, the default is 900x400
    show_plot if True: calling this function will show the plot, 
              if False: this function returns and object representing the plot.
                        it can then be fit into another bokeh layout object or shown with show()
    show_slider if True: the plot will feature a slider that when moved will hide datapoints after the shown date
                It is useful for looking at how good the models are at forecasting the next day
    '''
    
    if title is None: title = "Toxicity and Cell Count"
    p = figure(title = title, 
               x_axis_type = 'datetime', 
               x_axis_label = 'Date', 
               y_axis_label = 'ESP Estimate (cells/liter)',
               height = height, width = width)
    p.toolbar.active_scroll = p.select_one(WheelZoomTool)
    if hide_toolbar:
        p.toolbar.logo=None
        p.tools = []
    
    
    #Setting x_axis format
    if x_start and x_end: p.x_range = Range1d(start=x_start,end=x_end)
    elif x_start: p.x_range = Range1d(start=x_start)
    elif x_end:   p.x_range = Range1d(end=x_end)
    p.xaxis.formatter = DatetimeTickFormatter(days=['%Y-%m-%d'],months=['%Y-%m-%d'])

    # Setting axix ranges
    p.y_range = Range1d(start=0,end=4000)

    # Adding right-hand y axis
    p.extra_y_ranges = {"tox_axis": Range1d(end=300)}
    p.add_layout(LinearAxis(y_range_name="tox_axis", axis_label='Shellfish toxicity (µg/100g)'), 'right')

    # axis sizes
    p.title.text_font_size = '12pt'
    p.xaxis.axis_label_text_font_size  = '18pt'
    p.xaxis.major_label_text_font_size = '10pt'
    p.yaxis.axis_label_text_font_size  = '18pt'
    p.yaxis.major_label_text_font_size = '14pt'

    sources = dict()
    
    # plot each df column
    for series_name,c in zip(df,COLORS):

        no_nan_df = df[[series_name]].dropna(axis=0, how='any') 
        no_nan_xy = {'x':no_nan_df.index.tolist(), 'y':no_nan_df[series_name].tolist()}
        source = ColumnDataSource(no_nan_xy)
        sources[series_name] = source

        glyph_size = 6

        if 'esp' in series_name.lower() or 'ifcb ' in series_name.lower():     # cell counts, perhaps from different stations
            p.line('x','y', source = source, legend = series_name, line_color=c, line_width=glyph_size)        
            p.diamond('x','y', source = source, legend = series_name, color=c, line_color='black', size=glyph_size*2.2) 
        elif 'og copy MODEL' in series_name:
            p.circle('x','y', source = source, legend = series_name, 
                        y_range_name = 'tox_axis', color = c, size = 10)            
        elif 'og copy TOX' in series_name:
            p.circle('x','y', source = source, legend = series_name, 
                     y_range_name = 'tox_axis', color = 'white', line_color=c, size=10)              
        elif 'model' in series_name.lower() or '-' in series_name: # MODELed toxicity
            p.line('x','y', source = source, legend = series_name,
                   y_range_name = 'tox_axis', line_color = c, line_dash=[15,5], line_width=glyph_size)
            #p.circle('x','y', source = source, legend = series_name,
            #         y_range_name = 'tox_axis', color = c, line_color='black', size=glyph_size*2.2)

        else:                                 # observed/measured toxicity
            p.line('x','y', source = source, legend = series_name,
                   y_range_name = 'tox_axis', line_color = c, line_dash = 'solid', line_width=glyph_size)
            p.circle('x','y', source = source, legend = series_name,
                     y_range_name = 'tox_axis', color = c, line_color='black', size=glyph_size*2.2)
    
    #p.legend.location = "top_left"
    p.legend.click_policy="hide"
    
    # toxicity danger limit horizontal dotted line
    tox_danger_span = Span(location = 80, dimension = 'width',
                           line_color = 'grey',line_dash = 'dashed', line_width = 1,
                           y_range_name = 'tox_axis',tags=['Span','tox_danger_span'])
    p.add_layout(tox_danger_span)    

    ## SLIDER THINGS ##
    if show_slider:
        s,e = df.first_valid_index() , df.last_valid_index()
        s = pd.datetime(s.year,s.month,s.day)
        e = pd.datetime(e.year,e.month,e.day)
        slider = DateSlider(start=s,end=e,step=1,value=e,title="Date Slider", width= width-100, tags=['DateSlider','slider'])
        date_slider_span = Span(location=0,dimension='height', line_color='grey',line_dash='dashed', line_width=3, tags=['Span','date_slider_span'])
        p.add_layout(date_slider_span)

        update_slider = CustomJS(args=dict(sources=sources, span=date_slider_span, slider=slider), code="""
            span.location = slider.value

            for (var key in sources)
            {
                var data = sources[key].data
                var count = 0
                var selection = []
                for (var i in data['x'])
                {                                                // For each datapoint in each series (key)
                    selection.push(count)                        // add points to selection list
                    count = count + 1                            // until
                    if (data['x'][i] > slider.value) { break; }  // points go beyond slider-date value.
                }                                                // Remove the one extra selected point
                if ( key.toUpperCase().indexOf("MODEL") === -1 ) // -unless- it's the MODELed series'.                 
                    { selection.pop() }                          // Keep that one. 
                if (selection.length === 0)                      // An empty selection array (default) actually shows all points again.
                    { selection.push(0) }                        // so make sure to keep the first point selected to keep the rest hidden

                sources[key].selected.indices = selection        // Set the selection        
            }
            console.log(' ')
    """)
        slider.js_on_change('value', update_slider)
        layout = column(slider,p)
    else: layout = column(p)
    
    if show_plot:
        show(layout)
    else:
        return layout


#plot_shelltox(  df_test, show_plot=True, show_slider=True  )


# # The Model
#     T(t+1) = c1 * T(t) + c2 * C(t+1-tau)

# In[5]:


## THE MODEL ##
def model_tox_series(df, cell_col, tox_col, tox_col_model, tox_const, cell_const, lag, T0):
    '''
    Given a dataset and c1 c2 constants, generates modeled toxicity data 
    following the following equation for all known and interpolated cell_count values,
    Starting with the first measured toxicity value.
    T(t+1) = c1 * T(t) + c2 * C(t+1-tau)
    
    df is the data frame containing ESP and measured shellfish toxicity
    cell_col and tox_col are the column names of that ESP and shellfish toxicity data, respectively
    tox_col_model is the name to be given to the newly generated model data
    tox_const and cell_const are c1 and c2 respectivly. 
    c1 = 1-gamma, where gamma (γ) is the depuration constant.
    c2 = alpha, where alpha (α) is the uptake constant
    lag is the estimated time it takes for cells identified at an ESP to reach a downstream shellfish site,
        in integer days. Typically this value is 0 days or 1 days.
    T0 is the initial toxicity value. By default it is 0. 
        
    Returns: an rmsd value and the input dataframe including the tox_col_model data column.
    '''
    
    # Time Displacements
    aday = pd.Timedelta(days=1)
    lag = pd.Timedelta(days=lag)

    # Set up model column and interpolation
    if tox_col_model not in df.columns:
        df[tox_col_model] = pd.np.nan

    if tox_col == '':
        df_model = df[[cell_col,tox_col_model]].interpolate(limit_area='inside')
        if T0 is None: T0 = 0  # there can be no 'measured' start value
    else:
        df_model = df[[cell_col,tox_col,tox_col_model]].interpolate(limit_area='inside')

    # determining starting conditions
    for date in df_model.index:
        t0=date
        if T0 is None: # use measured value
            TO=df_model.loc[date,tox_col]
        else:          # use provided start value
            TO = T0
        C0=df_model.loc[date+aday-lag,cell_col]
        if not isnan(C0) and not isnan(TO): 
            df_model.loc[t0+aday,tox_col_model] = tox_const*TO + cell_const*C0
            break       

    # applying the model until cell counts run out
    for date in df_model[t0+aday:].index:
        Ti = df_model.loc[date,tox_col_model]
        Ci = df_model.loc[date+aday-lag,cell_col]
        if isnan(Ci): break
        df_model.loc[date+aday,tox_col_model] = tox_const*Ti + cell_const*Ci

    # calculating the root mean square deviation
    if tox_col != '':
        rms = calc_rms_df(df_model, tox_col, tox_col_model)
    else: rms = float('inf')

    # inserting the model data back into the input dataframe
    df[tox_col_model] = df_model[tox_col_model]
    
    return rms, df


# #### Root Mean Square Deviation
#     rmsd = √[ ∑(x-y)² / n ]
# where n is the total number of data points for x and y <br/>
# and where x and y are measured and modeled toxicity respectively, for a given date <br/>
# 
# rmsd is used over a simple sum of squares so as to give a basis of comparison between datasets of different lengths.
# 

# In[6]:


## Regression Line Functions ##

def calc_rms(x,y):
    '''
    Calculate the root mean square deviation for two lists of numbers, x and y.
    returns rms = √[ ∑(x-y)² / n ]
    '''
    xy = list(zip(*[x,y]))
    if len(xy) == 0: return float('inf')
    rmsd = sqrt(sum([(x-y)**2 for x,y in xy])/len(xy))
    return rmsd
    #squares = sum([(x-y)**2 for x,y in xy])
    #return squares

def calc_rms_df(df,measured_tox_col, modeled_tox_col):
    '''
    Calculate the root mean square deviation for two columns of a DataFrame
    returns √[ ∑(x-y)² / n ]
    '''
    df_rms = df[[measured_tox_col, modeled_tox_col]].dropna()
    toxobs = df_rms[measured_tox_col].tolist()
    toxmod = df_rms[modeled_tox_col].tolist()
    rms = calc_rms(toxobs,toxmod)
    return rms

def calc_rms_string(csv_string,measured_tox_col, modeled_tox_col):
    og_copy_df = pd.read_csv(StringIO(csv_string), parse_dates=True, index_col=0)
    rms = calc_rms_df(og_copy_df, measured_tox_col, modeled_tox_col)
    return rms 

def calc_regression(x,y):
    """
    y = m*x + b
    m = ∑(x-x̄)(y-ȳ)/∑(x-x̄)²
    b = ȳ - m*x̄
    
    Takes in two lists of numbers
    Returns a dict with 'slope', 'y_intercept', and 'r2'
    """
    xy = list(zip(*[x,y]))  # list of tuples
    x_bar = mean(x)
    y_bar = mean(y)
    xy_bar = mean([X*Y for X,Y in xy])
    x2_bar = mean([X**2 for X in x])
    y2_bar = mean([Y**2 for Y in y])
    
    rise = sum([(X-x_bar)*(Y-y_bar) for X,Y in xy])
    run  = sum([(X-x_bar)**2 for X,Y in xy])

    slope =  rise / run
    y_intercept = y_bar - slope*x_bar
    r2 = ( (xy_bar-x_bar*y_bar) / sqrt((x2_bar-x_bar**2)*(y2_bar-y_bar**2)) )**2
    rms = calc_rms(x,y)
    
    return {'slope':slope, 'y_intercept':y_intercept, 'r2':r2, 'rms':rms} 

def calc_r2_df(df,toxobs,toxmod):
    df_model = df[[toxobs, toxmod]].dropna()
    toxobs = df_model[toxobs].tolist()
    toxmod = df_model[toxmod].tolist()
    params = calc_regression(toxobs,toxmod)
    return params['r2']


def regression_plot(x,y, regression_params=None, show_plot=False):
    """
    Makes a scatterplot of x,y data and adds a trendline
    calc_regression() is used to calculate the trendline, 
    or the values can be passed in as a regression_params dict with 'slope', 'y_intercept', and 'r2'
    Returns the figure object p or if show_plot is True, displays the figure directly
    """
    
    # if regression_params not yet defined, calculate them
    if regression_params is None:
        regression_params = calc_regression(x,y)

    # extract regression params for easier access
    m = regression_params['slope']
    b = regression_params['y_intercept']
    r2 = regression_params['r2']
    rms = regression_params['rms']
    
    # Make the figure and plot the data
    p = figure(title='Toxicity Model vs. Observed',
            x_axis_label='Observed Toxicity', 
            y_axis_label='Modeled Toxicity',
            height = 250, width = 600)
    p.circle(x,y)

    # Add the trendline
    ## Bokeh version 13.0+ only ##
    slope = Slope(gradient=m, y_intercept=b,
                  line_color='orange', line_dash='dashed', line_width=2)
    p.add_layout(slope)
    
    # add trendline info on the plot, optionally.

    desc_text = 'Regression Line: y=mx+b for m={slope:.3f} b={y_intercept:.3f}  (r²={r2:.3f})'
    desc_text = desc_text.format(**regression_params)
    #p.text(x=[0],y=[0],text=[text])
    
    # show or return the figure
    if show_plot:
        print(desc_text)
        show(p)
    else:
        return desc_text,p,regression_params['r2']

def plot_many_regressions(df, measured_tox_col, modeled_tox_col):
    """
    Input: df, column name of measured toxicity series, column name of modeled toxicity series
    Returns 1) descriptive text of the regression line parameters
            2) bokeh Tabs layout of e regression line plot and data point DataTable
    """
    # setting up the data, keep only rows with both measured and modeled datapoints
    df_model = df[[measured_tox_col, modeled_tox_col]].dropna()
    toxobs = df_model[measured_tox_col].tolist()
    toxmod = df_model[modeled_tox_col].tolist()
    
    # DataTable tab 
    columns = [TableColumn(field='date', title='Date', formatter=DateFormatter(format="%Y-%m-%d")),
               TableColumn(field=measured_tox_col, title=measured_tox_col, formatter=NumberFormatter(format='0.0')),
               TableColumn(field=modeled_tox_col,  title=measured_tox_col, formatter=NumberFormatter(format='0.0'))]

    t = DataTable( columns=columns, source=ColumnDataSource(df_model), height=250)
    table_tab = Panel(child=t, title='Model Table')

    # Regression line tab
    desc_text, p, r2 = regression_plot(toxobs,toxmod)
    regline_tab = Panel(child=p, title='Model Trendline')

    # laying out the tabs
    tabbed_layout_output = Tabs(tabs=[regline_tab,table_tab], width=400)  
    
    return desc_text, tabbed_layout_output, r2
    
  
    


# ## Determining c1 c2
# A few methods were used to determine the best c1,c2 values for a given model's ESP and shellfish site. All methods involve iterating over a number c1,c2 values, calculating the root mean square deviation (rmsd) between the site's observed toxicity and the model, and keeping the c1,c2 values that produced the smallest rmsd.
# By default, the range of explored c1 and c2 values are 0.01-1 and 0.001-0.2 at a resolution of 0.01 and 0.001, respectively.
# The most appropriate lag (time of drift from ESP to shellfish toxicity site, in integer days, typically 0 or 1), is not determined by these functions and must be provided to calculate with. lag=0 is the default.
# 
# #### Brute Force
# The brute force approach is slowest but most thorough. It churns out 20,000 rmsd results (for the default ranges). Although slow, it's useful for saving the results as a csv for c1c2_heatmap() to plot the spread of good c1 c2 values over the whole range. For this function, an integer range CAN be supplied for lag.
# 
# #### 1D Recursive
# This method, dubbed semi-recursive method is faster than brute force thanks to its 1D binary search, but it's still too slow. For every c2 value it searches for the best c1 value using a binary search. 
# It could/should be upgraded to a full 2D binary min-value search algorithm which should prove to be very fast.
# 
# #### Recursive Walk
# The recursive walk method is the fastest and seems reliable enough to use as the primary method of c1,c2 calculating for interactive plots. This method recursively "walks downhill", calculating rmse values at and adjacent to a given c1,c2 location, picks the direction with the smallest resultant c1,c2 value, and repeats these steps with the new c1,c2 location. Each iteration brings the algorithm closer to the c1,c2 location with the minimum rmsd value. Once it finds a c1,c2 location with no lower adjacent rmsd values, it stops and reports the final c1,c2 location; a trough or reverse-peak in the c1,c2,rmsd map.
# By picking a starting location near where most of the best c1,c2 values lay (high/max c1,low/min c2 values), this algorithm so far proves to be fastest. This is in part due to the smoothly sloping nature of the c1,c2,rmsd map, which is better visualized later using the c1c2_heatmap() function. There you will even be able to see the path taken by the recursive walk function.
# 

# In[7]:


#Brute Force Method 

def model_autotune_brute(df, cell_col, tox_col, tox_model_col, tox_const_range=None, cell_const_range=None, lag_range=[0,1], T0=None, csv_filename=None):
    # calculate rms for all values in range of cell and tox constants
    if tox_const_range is None:
        tox_const_range = [round(n,5) for n in pd.np.linspace(1/100,1,100)]
    if cell_const_range is None:
        cell_const_range = [round(n,5) for n in pd.np.linspace(1/1000,200/1000,200)]
    
    # progress utility
    count = 0
    simu_len = len(cell_const_range)*len(tox_const_range)
    print('BRUTE FORCE COUNT:',simu_len)
    print('0.0%', end=' ')
    
    brute_force_results = []
    for lag in lag_range:
        for tox_const in tox_const_range:
            for cell_const in cell_const_range:
                rms,_ = model_tox_series(df, cell_col, tox_col, tox_model_col, tox_const, cell_const, lag, T0)
                brute_force_results.append((rms,tox_const,cell_const))
                
                #progress bar utility
                count +=1
                if count%50 ==0:       # Every 50 result, 
                    print('.',end='')  # show a progress dot
                if count%5000 == 0:    # Every 5000 results, give a completion percentage 
                    print(); print('{}%'.format(100*count/simu_len),end=' ')

    print()
    
    rms,tox_const,cell_const = min(brute_force_results)
    rms,best_df = model_tox_series(df, cell_col, tox_col, tox_model_col, tox_const, cell_const, lag, T0)

    # writing to csv, if csv_filename was given
    if csv_filename:
        pd.DataFrame(brute_force_results).to_csv(csv_filename, header=['rms','c1','c2'], index=False)
        
    return rms,tox_const,cell_const,best_df,brute_force_results


# In[8]:


# 1D recursive method

def model_autotune_semi_recursive(df, cell_col, tox_col, tox_model_col, tox_const_range=None, cell_const_range=None, lag=0, T0=0):
    # calculates rms for c1,c2 values using a binary search algorithm on c1 values for all possible c2 values
    if tox_const_range is None:
        tox_const_range = [round(n,5) for n in pd.np.linspace(1/100,1,100)]
    if cell_const_range is None:
        cell_const_range = [round(n,5) for n in pd.np.linspace(1/1000,200/1000,200)] 
    
    # The 1D recursive method
    def model_autotune_semi_recursive_tox(df, cell_col, tox_col, tox_model_col, cell_const, tox_const_range, lag, T0):

        # determin MID value for a given cell_const
        tox_mid_val_index = int(len(tox_const_range)/2)
        tox_mid_val = tox_const_range[tox_mid_val_index]
        rms_mid_val,df_mid = model_tox_series(df, cell_col, tox_col, tox_model_col, tox_mid_val,  cell_const, lag, T0)

        # determin MID+1 value for a given cell_const
        try:
            tox_next_val= tox_const_range[tox_mid_val_index+1]
            rms_next_val,_ = model_tox_series(df, cell_col, tox_col, tox_model_col, tox_next_val, cell_const, lag, T0)
        except IndexError: 
            rms_next_val= float('inf')

        # determine MID-1 value for a given cell_const
        try:
            if tox_mid_val_index-1 < 0: raise IndexError
            tox_prev_val= tox_const_range[tox_mid_val_index-1]
            rms_prev_val,_ = model_tox_series(df, cell_col, tox_col, tox_model_col, tox_prev_val, cell_const, lag, T0)    
        except IndexError: 
            rms_prev_val= float('inf')    


        # Minimum Value base case
        if rms_prev_val > rms_mid_val < rms_next_val:
            #congratulations, this is a local minimum!
            return rms_mid_val, tox_mid_val, cell_const, df_mid
        
        # Better c1 value is higher up, truncate away all lower c1 values and recurse
        # we want to look forwards to the later half of tox values for the rms minimum
        elif rms_prev_val > rms_mid_val > rms_next_val: 
            return model_autotune_semi_recursive_tox(df, cell_col, tox_col, tox_model_col, 
                                            tox_const_range=tox_const_range[tox_mid_val_index+1:], #truncation here
                                            cell_const=cell_const, lag=lag, T0=T0)
        # Better c1 value is lower down, truncate away all higher c1 values and recurse
        # we want to look backwards to the smaller half of tox values for the rms minimum        
        else:   
            return model_autotune_semi_recursive_tox(df, cell_col, tox_col, tox_model_col, 
                                            tox_const_range=tox_const_range[:tox_mid_val_index], #truncation here
                                            cell_const=cell_const, lag=lag, T0=T0)
    ## End of recursive function definition ##
    
    # Iteration over all c2 values, finding best c1 values with recursive binary search defined above.
    count=0
    semi_brute_results = []
    for cell_const in cell_const_range:
        rms,c1,c2,_ = model_autotune_semi_recursive_tox(df, cell_col, tox_col, tox_model_col, cell_const, tox_const_range, lag=lag, T0=T0)
        semi_brute_results.append((rms,c1,c2))
        
        # TODO, if the most recent rms is is bigger than the previous rms, stop right there! 
        #       it's probably the droid we're looking for. probably
        
        # progress bar printing
        count += 1
        print('.',end='')   # print a dot for every c2 row where a best c1 value was found
        if count%100 == 0:  # every 100 c2 rows, print the percentage of remainting values.
            print(); print('{}%'.format(100*count/len(cell_const_range)),end=' ')    
    
    # find the set with the smallest rms, and create the final best-fit model with it.
    rms,tox_const,cell_const = min(semi_brute_results)
    rms,best_df = model_tox_series(df, cell_col, tox_col, tox_model_col, tox_const, cell_const, lag, T0)
    
    return rms,tox_const,cell_const,best_df
        
    


# In[9]:


# Recursive Walk method - Fastest. must include diagonals to be accurate.

def model_autotune_recursive_walk(df, cell_col, tox_col, tox_model_col, path=None, lag=0, diag=True, T0=0): 
    
    # copy df so as not to affect copy outside of the recursive function. yes it's neccessary
    # todo: figure out how to avoid this to save on memory, cos seriously >_<
    df = df.copy()
    
    # we can afford to iterate of the full 0-1 span for c2 with this function
    c1_range = [round(n,5) for n in pd.np.linspace(1/100,1,100)]
    c2_range = [round(n,5) for n in pd.np.linspace(1/1000,1,1000)]

    ## NON RECURSIVE, INITIAL FUNCTION CALL CODE ##
    # Setting initial conditions and starting the recursion
    if path is None:
        path=[]
        # The following corner start locations tend to be closest to the peak
        ic1 = len(c1_range)-1 
        ic2 = 0 
        c1,c2 = c1_range[ic1], c2_range[ic2]
        rms_start,_= model_tox_series(df, cell_col, tox_col, tox_model_col, c1, c2, lag, T0)
        path.append((rms_start,c1,c2,ic1,ic2,'C'))
        rms_end, c1, c2, best_df, path = model_autotune_recursive_walk(df, cell_col, tox_col, tox_model_col, path, lag, diag, T0)
        return (rms_end, c1, c2, best_df, path)
    
    
    ### EVERYTHING BELOW HAPPENS DURING RECURSION ###
    
    # Taking stock of the algorithms current location in the c1,c2 map
    # path[-1] is the most recently added c1,c2 point during the previous iteration
    rms,c1,c2,ic1,ic2,_ = path[-1]
    
    # this lists the indexes of all the previously visited c1,c2 locations
    # we note them so as not to re-calculate any locations we've already been to
    ic12_path = [(p[3],p[4]) for p in path]

    
    ## Calculating Adjacent rms values ##
    
    # North, increasing ic1
    if ic1+1 < len(c1_range) and (ic1+1,ic2) not in ic12_path:
        c1N,c2N = c1_range[ic1+1], c2_range[ic2]
        rmsN,_= model_tox_series(df, cell_col, tox_col, tox_model_col, c1N, c2N, lag, T0)
    else: rmsN=float('inf')

    # East, increasing ic2
    if ic2+1 < len(c2_range) and (ic1,ic2+1) not in ic12_path:
        c1E,c2E = c1_range[ic1], c2_range[ic2+1]
        rmsE,_= model_tox_series(df, cell_col, tox_col, tox_model_col, c1E, c2E, lag, T0)
    else: rmsE=float('inf')

    # South, decreasing ic1
    if ic1-1 >= 0 and (ic1-1,ic2) not in ic12_path:
        c1S,c2S = c1_range[ic1-1], c2_range[ic2]
        rmsS,_= model_tox_series(df, cell_col, tox_col, tox_model_col, c1S, c2S, lag, T0)
    else: rmsS=float('inf')

    # West, decreasing ic2
    if ic2-1 >= 0 and (ic1,ic2-1) not in ic12_path:
        c1W,c2W = c1_range[ic1], c2_range[ic2-1]
        rmsW,_= model_tox_series(df, cell_col, tox_col, tox_model_col, c1W, c2W, lag, T0)
    else: rmsW=float('inf')
        
    ## Aggregating rms values ##
    directions = {'here':rms, 'north':rmsN, 'east':rmsE, 'south':rmsS, 'west':rmsW,
                  'ne':float('inf'), 'se':float('inf'), 'nw':float('inf'), 'sw':float('inf')}

    ## Calculating DIAGONALS ##
    # this helps the algorithm not get stuck in "mini-peaks"
    if diag:
        # North East diagonal
        if ic1+1 < len(c1_range) and ic2+1 < len(c2_range) and (ic1+1,ic2+1) not in ic12_path:
            c1NE,c2NE = c1_range[ic1+1], c2_range[ic2+1]
            rmsNE,_= model_tox_series(df, cell_col, tox_col, tox_model_col, c1NE, c2NE, lag, T0)
        else: rmsNE=float('inf')
        # South East diagonal
        if ic1-1 >= 0 and ic2+1 < len(c2_range) and (ic1-1,ic2+1) not in ic12_path:
            c1SE,c2SE = c1_range[ic1-1], c2_range[ic2+1]
            rmsSE,_= model_tox_series(df, cell_col, tox_col, tox_model_col, c1SE, c2SE, lag, T0)
        else: rmsSE=float('inf')
        # North West diagonal
        if ic1+1 < len(c1_range) and ic2-1 >= 0 and (ic1+1,ic2-1) not in ic12_path:
            c1NW,c2NW = c1_range[ic1+1], c2_range[ic2-1]
            rmsNW,_= model_tox_series(df, cell_col, tox_col, tox_model_col, c1NW, c2NW, lag, T0)
        else: rmsNW=float('inf')
        # South West Diagonal 
        if ic1-1 >= 0 and ic2-1 >= 0 and (ic1-1,ic2-1) not in ic12_path:
            c1SW,c2SW = c1_range[ic1-1], c2_range[ic2-1]
            rmsSW,_= model_tox_series(df, cell_col, tox_col, tox_model_col, c1SW, c2SW, lag, T0)
        else: rmsSW=float('inf')        
        directions.update({'ne':rmsNE, 'se':rmsSE, 'nw':rmsNW, 'sw':rmsSW})

    
    ## DETERMINE DIRECTION OF SMALLEST RMS ##
    downhill = min(directions, key=lambda key:directions[key])

    # Add the point in the determined direction to path and continue
    # TODO relabel the arrows to they match the directions taken on the heatmap. low priority.
    if downhill == 'north':
        print('↑',end='')
        path.append( (rmsN,c1N,c2N,ic1+1,ic2,'N') )
    elif downhill == 'east':
        print('→',end='')
        path.append( (rmsE,c1E,c2E,ic1,ic2+1,'E') )            
    elif downhill == 'south':
        print('↓',end='')
        path.append( (rmsS,c1S,c2S,ic1-1,ic2,'S') )
    elif downhill == 'west':
        print('←',end='')
        path.append( (rmsW,c1W,c2W,ic1,ic2-1,'W') )

    #Diagonals
    elif downhill == 'ne':
        print('↗',end='')
        path.append( (rmsNE,c1NE,c2NE,ic1+1,ic2+1,'ne') )                 
    elif downhill == 'se':
        print('↘',end='')
        path.append( (rmsSE,c1SE,c2SE,ic1-1,ic2+1,'se') )
    elif downhill == 'nw':
        print('↖',end='')
        path.append( (rmsNW,c1NW,c2NW,ic1+1,ic2-1,'nw') )       
    elif downhill == 'sw':
        print('↙',end='')
        path.append( (rmsSW,c1SW,c2SW,ic1-1,ic2-1,'sw') )

    ## THE BASE CASE! ##
    # generate the final model and return all the way up the recursion tree
    else: #downhill == 'here'
        print('⭑')
        rms,best_df = model_tox_series(df, cell_col, tox_col, tox_model_col, c1, c2, lag, T0)
        return (rms, c1, c2, best_df, path)
    
    # if the base case didn't trigger, re
    # repeats using the most recently added path point as the new starting point
    return model_autotune_recursive_walk(df, cell_col, tox_col, tox_model_col, path, lag, diag, T0)
        

# BRUTE FORCE MODEL Generation
# Change this from "Raw NBConvert" to "Code" to run (this may take like an hour to run them all)

# These models are ones we know the historic c1,c2 values used
#df2014 =load_dataset('2014.csv')
#
#print('2014 Jake')
#jake_result = model_autotune_brute(df2014, 'ESP Jake', 'TOX Potts Point', 'Jake-Potts', lag_range=[0], 
#                              csv_filename = 'c1c2maps/c1c2map.2014.Jake-Potts.csv')
#
#print('2014 Don')
#don_result = model_autotune_brute(df2014, 'ESP Don', 'TOX York River', 'Don-York', lag_range=[0],
#                                  csv_filename = 'c1c2maps/c1c2map.2014.Don-York.csv')
#
#print('2014 Dennis')
#dennis_result = model_autotune_brute(df2014, 'ESP Dennis', 'TOX Christmas Cove', 'Dennis-Christmas', lag_range=[1],
#                                     csv_filename = 'c1c2maps/c1c2map.2014.Dennis-Christmas.csv')
#
#print('2015 Dennis')
#df2015 = load_dataset('2015.csv')
#dennis15_result = model_autotune_brute(df2015, 'ESP Dennis', 'TOX Christmas Cove','Dennis-Christmas', lag_range=[1],
#                                       csv_filename = 'c1c2maps/c1c2map.2015.Dennis-Christmas.csv')
#
#print('2017 Roman')
#df2017 = load_dataset('2017.csv')
#roman_result = model_autotune_brute(df2017, 'ESP3 Roman', 'TOX Cutler Harbor', 'Roman-Culter', lag_range=[1],
#                                    csv_filename = 'c1c2maps/c1c2map.2017.Roman-Cutler.csv')
#
#print('DONE')

# # Comparing against old model data
# THE FOLLOWING FUNCTIONS ARE USED TO VALIDATE THE NEW MODEL by comparing it against the data created by the old model<br/>
# Note: og stands for original, as in the old c1,c2 values or the old model recorded in the spreadsheets

# In[10]:



def calculate_models(df, esp, shellfish_site, model_name, extra=None, og_c1=None,og_c2=None, c1=None,c2=None, lag=0, show_plot=True, T0=0): 
    """ 
    Generates and plots models using the original recorded og_c1,og_c2 values for comparison against auto-generated c1,c2 values
    extra, if provided, is a list of dates, original observed values, and original modeled values.
           These are then plotted along with the models generated here for comparison as large dots.
    
    Returns a dataframe with the generated models, and a plot object
            If show_plot is set to True, it only shows the plot onscreen, returning None.
    """
      
    # find best c1,c2 for use with new model if none were provided
    if c1 is None and c2 is None:
        #rms,c1,c2,df = model_autotune_semi_recursive(df, esp, shellfish_site, model_name,lag=lag,T0=T0)
        rms,c1,c2,df,path = model_autotune_recursive_walk(df, esp, shellfish_site, model_name,lag=lag, T0=T0)
    else: # or use pre-existing values
        rms,df = model_tox_series(df, esp, shellfish_site, model_name, c1, c2, lag=lag, T0=T0)
    print('        {} c1={} c2={} rmsd={:.3f}'.format(model_name.upper(),c1,c2,rms))

    # generate a model using old/historic values taken from spreadsheet
    if og_c1 is not None and og_c2 is not None:
        og_model_name = 'og c1c2 '+ model_name
        og_rms,df = model_tox_series(df, esp, shellfish_site, og_model_name, og_c1, og_c2, lag=lag, T0=T0)
        print('{} c1={:.2f} c2={:.3f} rmsd={:.3f}'.format(og_model_name,og_c1,og_c2,og_rms))
    else: og_model_name=''
    
    # columns to be returned at end of function
    col_list = [esp,shellfish_site,model_name,og_model_name]    
    
    # formats and merges fed-in "extra" data into the dataframe with the other series.
    if isinstance(extra,str):
        extra_df = pd.read_csv(StringIO(extra), parse_dates=True, index_col=0)
    elif isinstance(extra,list):
        extra_df = pd.DataFrame(extra[1:], columns = extra[0])
        extra_df['date'] = pd.to_datetime(extra_df['date'])
        extra_df.set_index('date',inplace=True)
    else:
        extra_df = extra
    if extra is not None:
        print('og copy {} c1=??   c2=??    rmsd={:.3f}'.format(model_name.lower(),
                                                               calc_rms_df(extra_df,
                                                                           extra_df.columns[1],
                                                                           extra_df.columns[0])))
        if all(col in df.columns for col in extra_df.columns):
            df[extra_df.columns] = extra_df[extra_df.columns]
        else:
            df = df.merge(extra_df, how='outer', on='date')
        col_list.extend(list(extra_df.columns))
    
    # plot the data
    p = plot_shelltox( df[col_list], title="Toxicity Model Comparisons" )
    
    # possible show plot and return data and plot object
    if show_plot:
        show(p)
    return df,p


def model_compare(df_tox, show_plot=False):
    """
    This function creates a tabbed interface showing 
    1)  a table with all the different models, but just on dates with data for all columns
    2+) trendline plots comparing the different models against the observed values
    """

    # setting up the math
    df_tox = df_tox.dropna()
    measured_tox_col = [col for col in df_tox.columns if col.upper().startswith('TOX') or col.upper().startswith('PSP')][0]
    modeled_toxs_cols = [col for col in df_tox.columns if not col.upper().startswith('TOX') and not col.upper().startswith('PSP')] 
    
    # setting up datatable columns
    measured_tox = df_tox[measured_tox_col].tolist()
    columns = [TableColumn(field='date', title='Date', formatter=DateFormatter(format="%Y-%m-%d")),
               TableColumn(field=measured_tox_col, title=measured_tox_col)]    
    
    # For each model, generate a tab with a trendline plot
    tabs = []
    for col in modeled_toxs_cols:
        modeled_tox = df_tox[col].tolist()
        columns.append(TableColumn(field=col, title=col))
        rms_regression_params = calc_regression(measured_tox, modeled_tox)
        desc_text,p,r2 = regression_plot(measured_tox, modeled_tox,rms_regression_params)
        text = PreText(text=desc_text, width=800)
        panel_child = column(p,text)
        regline_tab = Panel(child=panel_child, title=col+' Trendline')
        #regline_tab = Panel(child=p, title=col+' Trendline')
        tabs.append(regline_tab)

    # Generate and add DataTable tab 
    dt = DataTable( columns=columns, source=ColumnDataSource(df_tox), height=250)
    table_tab = Panel(child=dt, title='Model Table')
    tabs.insert(0,table_tab)

    # laying out the tabs
    tabbed_layout_output = Tabs(tabs=tabs, width=800)  
    
    # directly show or return the tabs 
    if show_plot:
        show(tabbed_layout_output)
    else:
        return tabbed_layout_output  


# In[11]:
# ~REDACTED~

# # Interactive Plots
# 
# While the plots generated above are interactive in the sense that a user can pan and zoom and move sliders and click legend items to affect the displayed data, all that functionality is merely baked into a static html document which includes clever javascript. 
# <br/><br/>
# To do more, like load specific datasets from a server or do calculations / generate data outside off the browser, one must set up a Bokeh server and write the interactive plot and features as a Bokeh "App". <br/>
# Fortunately, Bokeh and Jupyter Notebooks know how to work together on this.
# <br/><br/>
# With this functionality, it is straightforward to construct/prototype bokeh apps in a Notebook which can then be embedded on a website. This applies to interactive plots as well as the static plots we've previously seen.
# 
# ### multi_plot()
# multi_plot allows a user to select different csv datasets saved on the filesystem. Series of that dataset are all displayed by default, but the different series can be hidden or shown again by clicking on relevant series the figure's legend.
# 
# ### model_plot()
# model_plot is an inteactive model generator. First, select a dataset, and ESP series, and a Shellfish Site. Then either type in the desired c1,c2,lag values, or simply set the lag and let the program auto-calculate the best c1,c2 values for the given ESP and Shellfish Site! The model generated with those values will be automatically calculated and shown on the figure along with the ESP cell_counts and measured_toxicity values.<br/>
# A table and trendline will also be shown below the main figure to assess the quality of the model for the given sites and c1,c2 values. 

# In[12]:


##  MULTI_PLOT  ##

# Unique Imports #
from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application
from bokeh.models import Select
from bokeh.layouts import WidgetBox

import os

def modify_multi_plot(doc):
   
    # Load Initial Data
    global df
    datasets = sorted(os.listdir('data'))
    try: type(df)
    except NameError: df = load_dataset(datasets[-1])
        
    all_series = df.columns.tolist()
    esp_series = [serie for serie in all_series if 'esp' in serie.lower() or 'ifcb' in serie.lower()]
    tox_series = [serie for serie in all_series if serie.upper().startswith('TOX') or serie.upper().startswith('PSP') ]
    
    # Controls #
    year_select = Select(title='Dataset', options=datasets, value=datasets[-1])
    
    # Control Callbacks #
    def update_dataset(attr, old, new):
        global df
        df = load_dataset(new)
        fig = plot_shelltox(df)
        layout.children[1] = fig
        
        # changes to series_select.value triggers update_series which in turn updates the plot
        all_series = df.columns.tolist()
        series_select.value = []
        series_select.options = all_series
    year_select.on_change('value',update_dataset)
    
    # PUT PLOT & LAYOUT TOGETHER #
    fig = plot_shelltox(df)
    layout = column(year_select,fig)
    doc.add_root(layout)

def multi_plot():
   
    # Set up the Application 
    handler = FunctionHandler(modify_multi_plot)
    app = Application(handler)

    # Create the Document
    # Not strictly necessary, but helps w/ debugging
    doc = app.create_document()

    # Show the application
    # Make sure the URL matches your Jupyter instance
    show(app, notebook_url="localhost:8888")

#multi_plot()


# In[13]:


## INTERACTIVE MODEL EXPLORER ##

# Unique Imports #
from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application
from bokeh.models import Select, TextInput, PreText, Button
import os

def modify_model_plot(doc):
    '''This function is the heart of an interactive online "app" plotting tool'''
    T0 = None # initial toxicity of model will be first available measured toxicity
    T0 = 0    # initial toxicity of model will be zero
    
    # Load Initial Data
    global df
    datasets = sorted(os.listdir('data'))
    try: 
        df = load_dataset(datasets[-1])
        load_error = False 
    except: 
        df = pd.DataFrame()
        load_error = True
        
    all_series = df.columns.tolist()
    esp_series = [serie for serie in all_series if serie.upper().startswith('ESP') or serie.upper().startswith('IFCB')]
    tox_series = [serie for serie in all_series if serie.upper().startswith('TOX') or serie.upper().startswith('PSP') ]     
    
    # Control Widgets #
    year_select      = Select(title='Dataset',          options=datasets,   value=datasets[-1])
    model_esp_select = Select(title='ESP for Modeling', options=esp_series, value=esp_series[0] if esp_series else '')
    model_tox_select = Select(title='Shellfish Site',   options=tox_series, value=tox_series[0] if tox_series else '')
    cell_const_input = TextInput(title='Uptake Constant (α)',     value='0.014', placeholder='0.014')
    tox_const_input  = TextInput(title='Depuration Constant (γ)', value='0.1',   placeholder='0.1')
    lag_const_input  = TextInput(title='lag (integer days)',      value='0')
    autotune_button  = Button(label='Calculate α & γ', button_type="success")    
    autotune_alpha_button = Button(label='Calculate α (hold γ)', button_type="success")    
    export_button    = Button(label='Export to csv', button_type="success")
    tox_init_input   = TextInput(title='model initial toxicity (µg/100g)', value='', placeholder='most recent toxicity (interpolated)')
    
    # Create main plot and layout #
    fig = plot_shelltox(df, show_slider=False)
    layout = row(WidgetBox(year_select, model_esp_select, model_tox_select,
                           lag_const_input, tox_init_input,
                           cell_const_input, tox_const_input, 
                           autotune_button, 
                           autotune_alpha_button,
                           export_button,
                           width = 250),
                 column(fig, width = 800), 
                 height=850)
    
    # associate the features created above with online "app" document
    doc.add_root(layout)
    
    
    
    ### WIDGET CALLBACKS ###
    # these are allowed be defined after features have been added to the document root
    
    def update_dataset(attr, old, new):
        ''' update_dataset loads a csv file of date,ESP,TOX data from the ./data directory '''
        global df
        df = load_dataset(new)
        
        # when a new dataset is loaded, update the selector boxes trigger plot update on final update
        all_series = df.columns.tolist()
        esps = [serie for serie in all_series if serie.upper().startswith('ESP') or serie.upper().startswith('IFCB')]
        toxs = [serie for serie in all_series if serie.upper().startswith('TOX') or serie.upper().startswith('PSP')]
        model_esp_select.remove_on_change('value',update_plot) # prevents pre-mature triggering of
        model_tox_select.remove_on_change('value',update_plot) # tox and esp const update() callback
        model_esp_select.options = esps
        model_tox_select.options = toxs 
        model_esp_select.value = esps[0]  
        model_tox_select.value = toxs[0] if toxs else ''
        model_esp_select.on_change('value',update_plot)
        model_tox_select.on_change('value',update_plot)
        p = plot_shelltox(df, show_slider=False, title=new)
        layout.children[1] = column(p)
        export_button.callback = CustomJS(args=dict(file_content=df.to_csv(), file_name=new), 
                                          code=open(os.path.join(os.path.dirname(__file__), "download.js")).read())    
    year_select.on_change('value',update_dataset)

    def format_model_name(esp,tox):
        if 'IFCB' in esp.upper():
            esp_name = 'IFCB'
        else:
            esp_name = esp.split(' ')[1]
               
        try: tox_name = tox.split(' ')[1]
        except IndexError: tox_name=tox
        if tox_name == '':
            return '{} Model'.format(esp_name)
        else:
            model_name = '{}-{} Model'.format(esp_name, tox_name)
            return model_name
        

    def update_plot(attr, old, new):
        ''' update_plot creates a plot for the currently selected ESP and shellfish-site.
            If tox_const, cell_const, and lag are specified and valid, it adds a model to 
            the plot and show a trendline for the model.
        '''
        global df
        try: 
            tox_const = 1-float(tox_const_input.value) #c1
            cell_const = float(cell_const_input.value) #c2
            lag_const = int(lag_const_input.value)
            tox_init = tox_init_input.value
            if tox_init=='': tox_init=None # replaces T0 value
            else: tox_init = float(tox_init)

        except:
            plot_sans_model = plot_shelltox(df[[model_esp_select.value,model_tox_select.value]], show_slider=False)
            layout.children[1] = column(plot_sans_model)
            return

        model_name = format_model_name(model_esp_select.value,model_tox_select.value)

        # Calculate model points
        rms,df2plot = model_tox_series(df, model_esp_select.value, model_tox_select.value, model_name, 
                                       tox_const, cell_const, lag_const, tox_init)
        if model_tox_select.value != '':
            df2plot = df2plot[[model_esp_select.value,model_tox_select.value,model_name]]
        else:
            df2plot = df2plot[[model_esp_select.value,model_name]]

        # Calculate Model Stats
        try:
            descriptive_text,regression_tabbed_tableplot, r2 = plot_many_regressions(df2plot, model_tox_select.value, model_name)
            descriptive_text = PreText(text = descriptive_text, width=600)
        except: r2=0

        #Create Plot and update/replace the layout
        title = '{} at {} (α={:.3f}, γ={:.2f}, rmse={:.3f}, r2={:.1f}%)'.format(model_esp_select.value, 
                                                                                model_tox_select.value.replace('TOX ',''), 
                                                                                cell_const, 1-tox_const, rms, 100*r2)        
        p = plot_shelltox(df2plot, show_slider=False, title=title)
        try: layout.children[1] = column(p, regression_tabbed_tableplot, descriptive_text)
        except NameError: layout.children[1] = column(p)
        autotune_button.label = 'Calculate α & γ'
        autotune_alpha_button.label = 'Calculate α (hold γ)'
        
        export_filename = '{}_{}.csv'.format(year_select.value.replace('.csv',''), model_name.replace(' ','_'))
        export_button.callback = CustomJS(args=dict(file_content=df2plot.to_csv(), file_name=export_filename), 
                                          code=open(os.path.join(os.path.dirname(__file__), "download.js")).read())        
        
    # Widgets that on-change explicitely update the plot
    model_esp_select.on_change('value', update_plot)
    model_tox_select.on_change('value', update_plot)
    tox_const_input.on_change( 'value', update_plot)
    cell_const_input.on_change('value', update_plot)
    lag_const_input.on_change( 'value', update_plot)
    tox_init_input.on_change(  'value', update_plot)

    def autotune():
        ''' autotune calls a recursive function that determines the best c1 and c2 values
            for the selected ESP and shellfish-site.
        '''
        global df
        autotune_button.label = 'processing...'
        try: lag = int(lag_const_input.value)
        except ValueError: 
            print('LAG INVALID')
            autotune_button.label = 'error: lag invalid'
            return
        
        try:
            tox_init = tox_init_input.value.strip()
            if tox_init=='': tox_init=None # replaces T0 value
            else: tox_init = float(tox_init)
        except ValueError: 
            print('TOX_INIT INVALID') 
            autotune_button.label = 'error: initial_toxicity invalid'
            return
        
        model_name = format_model_name(model_esp_select.value, model_tox_select.value)
        rms, c1, c2, df, _ = model_autotune_recursive_walk(df, model_esp_select.value, 
                                                           model_tox_select.value, 
                                                           model_name, lag=lag, T0=tox_init)
                                                                   
        autotune_button.label = 'Calculate α & γ'
        tox_const_input.value, cell_const_input.value = str(round(1-c1,3)), str(round(c2,3))

    autotune_button.on_click(autotune)

    def autotune_alpha():
        ''' autotune calls a recursive function that determines the best c2 value while holding c1 (depuration) constant
            for the selected ESP and shellfish-site.
        '''
        global df
        autotune_alpha_button.label = 'processing...'
        try: lag = int(lag_const_input.value)
        except ValueError: 
            print('LAG INVALID')
            autotune_alpha_button.label = 'error: lag invalid'
            return
        
        try:
            tox_init = tox_init_input.value.strip()
            if tox_init=='': tox_init=None # replaces T0 value
            else: tox_init = float(tox_init)
        except ValueError: 
            print('TOX_INIT INVALID') 
            autotune_alpha_button.label = 'error: initial_toxicity invalid'
            return
        
        tox_const = 1-float(tox_const_input.value) #c1
        model_name = format_model_name(model_esp_select.value, model_tox_select.value)
        rms, c1, c2, df, _ = model_autotune_brute(df, model_esp_select.value, model_tox_select.value, model_name, 
                                                  lag_range=[lag], T0=tox_init, tox_const_range=[tox_const])

        autotune_alpha_button.label = 'Calculate α (hold γ)'
        tox_const_input.value, cell_const_input.value = str(round(1-c1,3)), str(round(c2,3))

    autotune_alpha_button.on_click(autotune_alpha)

    # initial download button callback. this callback is overwritten with the correct df every plot update
    export_button.callback = CustomJS(args=dict(file_content=df.to_csv(), file_name=datasets[-1]), 
                                      code=open(os.path.join(os.path.dirname(__file__), "download.js")).read())

    
def model_plot(notebook=True):
    # Set up the Application 
    handler = FunctionHandler(modify_model_plot)
    app = Application(handler)

    # Create the Document
    # Not strictly necessary, but helps w/ debugging
    doc = app.create_document()

    # Show the application
    # Make sure the URL matches your Jupyter instance
    if notebook:
        show(app, notebook_url="localhost:8888")      

#model_plot()



# # c1 c2 Heatmaps
# 
# By using data generated with the brute force c1,c2 determinating function we can get a sense of where and how rmsd is distributed in relation the the c1,c2 values picked to model for a given ESP and Shellfish site. <br/>
# This makes for an intuitive way to visualize where the best c1,c2 values reside specifically and generally. 

# In[14]:


## c1 c2 HEATMAP ##

from bokeh.palettes import Spectral11,Viridis11,Plasma256,Magma256,Viridis256,Inferno256
from bokeh.transform import linear_cmap

def c1c2_heatmap(fname, show_plot=True,high_res_color=True):
    '''
    fname is the file name of a csv with rms,c1,c2 columns
    show_plot if True directly shows the plot (default) 
              if False this function returns the plot object which can be modified and later shown
    high_res_color: if True  a 256 value color spectrum is used, creating a smooth blend of colors
                    if False, an 11 value color spectrum is used, creating a color banded, regionalized mapping
                    
    You can hover your mouse over the plot to see the c1,c2 values and their associated rmsd value. 
    You might have to zoom as when zoomed out the datapoint glyphs overlap.
    '''
    # read csv data
    df = pd.read_csv(fname)

    # find the row index containing the smallest rmsd value and extract the c1,c2,rmsd values from it
    imin = df['rms'].idxmin()
    min_rms,min_c1,min_c2 = df.iloc[imin].tolist()
    
    # formatting the plot title
    title = '{}: min(rms)={:.3f} at c1={:.2f} c2={:.3f}'
    title = title.format(fname,min_rms,min_c1,min_c2)
    
    # Defining the hover tooltip
    ttip = [("(c1,c2)", "(@c1{(.00)}, @c2{(.000)})"),
            ("rms", "@rms"),]
    
    # Creating the figure
    p = figure(title=title, tooltips = ttip)
    p.xaxis.axis_label = 'c1'
    p.yaxis.axis_label = 'c2'
    
    
    ## Color mapping ##
    # We want to better differentiate between small values, so lets take the inverse of rmsd values
    # otherwise it highlights the peaks instead of the troughs.
    df['1/rms'] = 1/df['rms']
    rms_inv = df['1/rms'].tolist()
    
    if high_res_color: pallette = Viridis256
    else:              pallette = Viridis11
        
    # linear_cmap is a color mapping function from bokeh.         
    mapper = linear_cmap(field_name='1/rms', palette= pallette,low=min(rms_inv) ,high=max(rms_inv))

    # Plotting all the c1,c2 values, color coded by rmsd mapper
    source = ColumnDataSource(df)    
    p.square(x='c1',y='c2',source=source,color=mapper, size=4)
    
    # highlighting/adding in red the smallest, most desirable rmsd value's c1,c2 location
    p.circle(x=min_c1,y=min_c2,color='red')
    
    #show or return plot
    if show_plot:
        show(p)
    else:
        return(p)


