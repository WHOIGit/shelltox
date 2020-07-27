#!/usr/bin/env python3
"""
A script for calculating the best uptake and depuration constants given historical data as input.
"""
import argparse
import pandas as pd
import numpy as np
from math import isnan,sqrt
from statistics import mean
import os


def calc_rms(x,y):
    """
    Calculate the root mean square deviation for two lists of numbers, x and y.
    returns rms = √[ ∑(x-y)² / n ]
    rms is used over a simple sum of squares so as to give a basis of comparison between datasets of different lengths.
    """
    xy = list(zip(*[x,y]))
    rms = sqrt(sum([(x-y)**2 for x,y in xy])/len(xy))
    return rms

def calc_sos(x,y):
    """
    Calculate the sum-of-squared error for two lists of numbers, x and y.
    returns sos = ∑(x-y)²
    """
    xy = list(zip(*[x,y]))
    sos = sum([(x-y)**2 for x,y in xy])
    return sos

def calc_rms_df(df, measured_tox_col, modeled_tox_col):
    """
    Calculate the root mean square deviation for two columns of a DataFrame
    returns √[ ∑(x-y)² / n ]
    """
    df_rms = df[[measured_tox_col, modeled_tox_col]].dropna()
    toxobs = df_rms[measured_tox_col].tolist()
    toxmod = df_rms[modeled_tox_col].tolist()
    rms = calc_rms(toxobs,toxmod)
    return rms
    
    

def load_dataset(fname,src=None):
    """Loads a csv with date, ESP, and TOX data columns"""
    
    # Ingest the csv, interpretting the first column as a dates index
    if src is not None:
        fname = os.path.join(src,fname)
    df = pd.read_csv(fname, index_col=0, parse_dates=True).dropna(how="all")

    # insert rows where there are missing days to get a monotonic timeseries
    df = df.reindex( pd.date_range(start=df.index[0],end=df.index[-1]) )
    df.index.name = 'date'
    return df



## THE MODEL ##
def model_tox_series(df, cell_col, tox_col, tox_col_model, tox_const, cell_const, T0=None, lag=0):
    """
    Given a dataset and c1 c2 constants, generates modeled toxicity data 
    following the following equation for all known and interpolated cell_count values,
    Starting with the first measured toxicity value.
    T(t+1) = c1 * T(t) + c2 * C(t+1-tau)
    
    df is the data frame containing ESP and TOX measured shellfish toxicity
    cell_col and tox_col are the column names of that ESP and shellfish toxicity data, respectively
    tox_col_model is the name to be given to the newly generated model data
    tox_const and cell_const are c1 and c2 respectivly. 
    tox_const = c1 = 1-gamma, where gamma is the depuration value
    cell_const = c2 = beta, where beta is the uptake value
    lag is the estimated time it takes for cells identified at an ESP to reach a downstream shellfish site,
        in integer days. Typically this value is 0 days or 1 days.
    T0 is the initial model toxicity value. 
        If None, the measured or interpolated toxicity value from 1 day prior to the first esp datum is used.
        
    Returns: an rms value, 
             the input dataframe with additional tox_col_model column, 
             list of measured tox data, 
             list of modeled tox data.
    """
    
    # Time Displacements
    aday = pd.Timedelta(days=1)
    lag = pd.Timedelta(days=lag)

    # Set up model column and interpolation
    df = df[[cell_col,tox_col]].interpolate(limit_area='inside')
    df[tox_col_model] = pd.np.nan
    
    # determining starting conditions
    for date in df.index:
        t0=date
        if T0 is None: # use measured value
            TO=df.loc[date,tox_col]
        else:          # use provided start value
            TO = T0
        C0=df.loc[date+aday-lag,cell_col]
        if not isnan(C0) and not isnan(TO): 
            df.loc[t0+aday,tox_col_model] = tox_const*TO + cell_const*C0
            break       
    
    # applying the model until cell counts run out
    for date in df[t0+aday:].index:
        Ti = df.loc[date,tox_col_model]
        Ci = df.loc[date+aday-lag,cell_col]
        if isnan(Ci): break
        df.loc[date+aday,tox_col_model] = tox_const*Ti + cell_const*Ci
    
    # calculating the root mean square deviation
    rms = calc_rms_df(df, tox_col, tox_col_model)
    
    # reporting list of real and modeled toxicity values for cumu_rms
    df_tox = df[[tox_col, tox_col_model]].dropna()
    tox_real = df_tox[tox_col].tolist()
    tox_model = df_tox[tox_col_model].tolist()
    #tox_tups = list(df[[tox_col, tox_col_model]].dropna().itertuples(index=False, name=None))

    return rms, df, tox_real, tox_model
    
    
    
def calculate_RMS_table(input_tups,beta_gammas, T0=None, output='df'):
    """
    Calculates RMS for all year-location-pairs for all beta_gammas. 
    Additionally calculates cumulative RMS across year-location-pairs for all beta_gammas.
    Input tups must be a list of tuples with 4 elements: model_column_label,esp_column_label,tox_column_label,dataframe
    beta_gammas is a list of beta and gamma tuples
    T0 is the initial model toxicity value. Default is None (select value automatically)
    output changes the output format. If 'df', a dataframe is returned. 
           If any other string, a csv with that filename is created.
           If anything else, eg None, a dict is returned.
    """
    print('Calculating RMS')
    betas,gammas = list(zip(*beta_gammas))
    df_data = dict(betas=betas,gammas=gammas)
    
    cumu_toxlist = dict()
    for a,g in beta_gammas:
        cumu_toxlist[(a,g)] = dict(real=[],model=[])
    
    for col,esp,tox,sub_df in input_tups:
        print('   ',col)
        beta_gamma_results = []
        for a,g in beta_gammas:
            rms,_,toxlist_real,toxlist_model = model_tox_series(sub_df, esp, tox, col, tox_const=1-g, cell_const=a, T0=T0)
            beta_gamma_results.append(rms)
            cumu_toxlist[(a,g)]['real'].extend(toxlist_real)
            cumu_toxlist[(a,g)]['model'].extend(toxlist_model)
        df_data[col] = beta_gamma_results

    print('   ','Cumulative RMS (cumu_rms)')
    for a,g in beta_gammas:
        cumu_toxlist[(a,g)]['rms'] = calc_rms(cumu_toxlist[(a,g)]['real'],cumu_toxlist[(a,g)]['model'])
        #cumu_toxlist[(a,g)]['sos'] = calc_sos(cumu_toxlist[(a,g)]['real'],cumu_toxlist[(a,g)]['model'])
    df_data['cumu_rms'] = [cumu_toxlist[(a,g)]['rms'] for a,g in beta_gammas]

    print()
    if output=='df':
        df = pd.DataFrame(df_data)
        return df.set_index(['betas','gammas'])
    elif isinstance(output,str):
        df = pd.DataFrame(df_data)
        df = df.set_index(['betas','gammas'])
        df.to_csv(output)
    else:
        return df_data


   
def annotate(df, outfile=None, print_final=[1,2,3]):
    """
    Adds footer and final columns with value means to input dataframe df
         Footer includes best-rms on a per-column basis and the associated beta and gamma values
    If outfile is specified, the table is saved to the entered filename.
       outfile may also be "stdout" in which case the a reduced summary table is displayed
    print_final is optional and (if outfile is not "stdout") shows the final results for method1, method2, and method3.
    Returns the modified df.
    """
    print('Annotating Table...')
    df = df.copy()
    #df = df.set_index(['betas','gammas'])
    cumu_series = df['cumu_rms']
    df = df.drop('cumu_rms',axis=1)
        
    valmin_mean1 = df.min().mean()
    beta_mins,gamma_mins = list(zip(*df.idxmin().tolist()))
    betamin_mean1,gammamin_mean1 = mean(beta_mins),mean(gamma_mins)

    df[''] = pd.Series()
    df['method1_mean'] = pd.Series()
    df['method2_mean'] = df.mean(axis=1)
    df['method3_cumu'] = cumu_series
        
    beta_mins,gamma_mins = list(zip(*[tup if isinstance(tup,tuple) else (np.nan,np.nan) for tup in df.idxmin().tolist()]))

    df = df.append(pd.Series(name=('','')))
    df.loc[('best','rms'),:] = df.min()
    df.loc[('best','beta'),:] = beta_mins
    df.loc[('best','gamma'),:] = gamma_mins

    df.loc[('best','rms'),'method1_mean'] = valmin_mean1
    df.loc[('best','beta'),'method1_mean'] = betamin_mean1
    df.loc[('best','gamma'),'method1_mean'] = gammamin_mean1

    valmin_mean2   = df.loc[('best','rms'),'method2_mean']
    betamin_mean2 = df.loc[('best','beta'),'method2_mean']
    gammamin_mean2 = df.loc[('best','gamma'),'method2_mean']
    
    valmin_cumu   = df.loc[('best','rms'),'method3_cumu']
    betamin_cumu = df.loc[('best','beta'),'method3_cumu']
    gammamin_cumu = df.loc[('best','gamma'),'method3_cumu']
    
    print('Outputting table to:', outfile)
    if outfile=='stdout':
        df_short = df.loc['best'].T.fillna('')
        df_short.columns.name=''
        print(df_short)
    elif outfile:
        df.to_csv(outfile)

    if print_final: print('Final Results:')
    if print_final==[3]:
        print('    least-rms={:.3f}, beta-gamma=({:.3f},{:.3f})'.format(valmin_cumu,betamin_cumu,gammamin_cumu))
    elif print_final:
        if 1 in print_final: print('  Method1:  mean-rms={:.3f}, beta-gamma=({:.3f},{:.3f})'.format(valmin_mean1,betamin_mean1,gammamin_mean1))
        if 2 in print_final: print('  Method2:  mean-rms={:.3f}, beta-gamma=({:.3f},{:.3f})'.format(valmin_mean2,betamin_mean2,gammamin_mean2))
        if 3 in print_final: print('  Method3: least-rms={:.3f}, beta-gamma=({:.3f},{:.3f})'.format(valmin_cumu,betamin_cumu,gammamin_cumu))

    return df        



def parse_input_args(targets, beta_arg, gamma_arg, src=None):
    """
    Converts cli input values to a format usable by calculate_RMS_table. 
    Also checks that input args are valid.
    targets is a list of strings where each string has the followwing format: "YEARfile,ESP column,TOX column"
    beta_arg and gamma_arg are either a single-item list containing a float, or a 3 item list with 3 floats: START STOP STEP
    """
    print('validating input...')
    input_tups = []
    for target in targets:
        print('   ',target)
        fname,esp,tox = target.split(',')
        if not fname.lower().endswith('.csv'): fname = fname+'.csv'
        if src is not None:
            fname = os.path.join(src,fname)
        year = os.path.splitext(os.path.basename(fname))[0]
        col = '{}:{}:{}'.format(year,esp,tox.replace('TOX ','')).replace(' ','_')
        assert os.path.isfile(fname), '{} was not found'.format(fname)
        sub_df = load_dataset(fname)
        assert esp in sub_df.columns, '"{}" is not a valid column in {}'.format(esp, fname)
        assert tox in sub_df.columns, '"{}" is not a valid column in {}'.format(tox, fname)
        sub_df = sub_df[[esp,tox]]
        input_tups.append((col,esp,tox,sub_df))
    
    print()
    if len(beta_arg)==3:
        start,stop,step = beta_arg
        betas = np.arange(start,stop+step*0.999,step)
        betas = [round(a, 6) for a in betas]
        print('Uptake Range: {} to {} ({} steps)'.format(start,stop,len(betas)))
    else:
        betas = beta_arg
        print('Uptake = {}'.format(betas[0]))
    
    if len(gamma_arg)==3:
        start,stop,step = gamma_arg
        gammas = np.arange(start,stop+step*0.999,step)
        gammas = [round(g, 6) for g in gammas]
        print('Depuration Range: {} to {} ({} steps)'.format(start,stop,len(gammas)))
    else:
        gammas = gamma_arg
        print('Depuration = {}'.format(gammas[0]))

    beta_gammas = []
    for a in betas:
        for g in gammas:
            beta_gammas.append((a,g))
    print()
    return input_tups,beta_gammas
    
               


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", nargs="+", metavar="YEAR,ESP,TOX", help='List of "YEAR,ESP,TOX" labels. "YEAR" is the filename of a YEAR.csv file found in SRC, "ESP" and "TOX" are column-headers from YEAR.csv. Note that each triplet should be space-delimited and surrounded by quotes "". This positional argument will also accept a configuration file in the described format, newline-delimited.')
    parser.add_argument("--src", help="Path to directory with input csv datafiles. Default is the current working directory.")

    groupB = parser.add_mutually_exclusive_group()
    groupB.add_argument("-b", "--betas", dest='beta_args', default=[0,0.05,0.001], 
                        nargs=3, type=float, metavar=('START','STOP','STEP'),
                        help='Range of uptake constants to asses. Default is "0 0.05 0.001". Mutually-exclusive with UPTAKE.')
    groupB.add_argument("-u","--uptake", metavar='UPTAKE',dest='beta_args', nargs=1, type=float, help='Single locked-in uptake constant. Mutually-exclusive with --betas.')
    
    groupG = parser.add_mutually_exclusive_group()
    groupG.add_argument("-d","--depuration",metavar='DEPURATION',dest='gamma_args',nargs=1,type=float, default=[0.1],
                        help='Single locked-in depuration constant. Default is "0.1". Mutually-exclusive with --gammas.')
    groupG.add_argument("-g", "--gammas", dest='gamma_args', #default=[0,0.2, 0.005], 
                        nargs=3, type=float, metavar=('START','STOP','STEP'),
                        help='Range of depuration constants to asses. Eg: "0 0.2 0.005". Mutually-exclusive with DEPURATION.')

    parser.add_argument("-o","--outfile", default='stdout', help='The file to output the final csv table to. If not specified, a summary table is displayed.')    

    parser.add_argument("-i","--model-init", metavar='I', type=float, default=None, 
                        help='If included, sets the initial modeled toxicity to the specified value, eg 0. If not included (default), the measured or interpolated toxicity from 1-day prior to the first ESP datapoint is used.')
    parser.add_argument("-m", "--method", metavar='N', type=int, nargs="+", choices=[1,2,3],
                        help=argparse.SUPPRESS)
                        #help='Determins the final shown command-line output text. The output csv includes all 3 methods of determining best uptake and depuration values, this merely determines what is shown onscreen when the process completes. For example "1 2 3" shows the final result using all 3 methods. If OUTFILE is specified, default is "3".')

    args = parser.parse_args()    

    if args.outfile!='stdout' and args.method is None:
        args.method=[3]

    if os.path.isfile(args.input[0]):
        config_file = args.input[0]
        print('loading from inputs from file: ' + config_file)
        args.input = []
        with open(config_file) as f:
            for line in f:
                line=line.strip()
                if line.replace(',','')=='': continue
                args.input.append(line.strip())

    ## Parse and Validate Inputs
    input_tups, beta_gammas = parse_input_args(args.input, args.beta_args, args.gamma_args, args.src)

    ## Calculate Results
    df = calculate_RMS_table(input_tups, beta_gammas, T0=args.model_init)
    df = annotate(df, outfile=args.outfile, print_final=args.method)
    
    
