from shelltox_api import modify_model_plot
import pull_new_esp_data
from bokeh.plotting import curdoc



print('START')

try: 
    print('updating csv from science.whoi.edu/fieldcelldata')
    pull_new_esp_data.do_default()
except Exception as e:
    print(type(e),e) 

modify_model_plot(curdoc())

print('ok loaded !!!')
