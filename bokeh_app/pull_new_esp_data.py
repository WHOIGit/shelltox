#!/usr/bin/env python
import pandas as pd
import requests
from io import StringIO
import argparse,os
from datetime import datetime

DEFAULT_ESPs = ['Dennis','Don','Jake','Roman']

def get_web_df(year, esps):
    url_template = "http://science.whoi.edu/fieldcelldata/secure/{year}/{esp}/counts.csv"
    user,pword = 'ESP_lab','ESPdata_GOM'
    cell_col = 'ESPvalue_cellsLiter'

    web_df = pd.DataFrame(pd.DataFrame(columns=[], index=pd.to_datetime([])))
    web_df.index.name = 'Date'

    for esp in esps:
        url = url_template.format(year=year,esp=esp.lower())
        resp = requests.get(url, auth=(user, pword))
        if resp.status_code != 200:
            raise requests.HTTPError(resp.status_code,url)
        csv_fio = StringIO(resp.text)
        df = pd.read_csv(csv_fio)
        df['Date'] = pd.DatetimeIndex(df.date).normalize()
        df = df.set_index('Date')
        
        date_dupes = df.index.duplicated(keep=False)
        #if any(date_dupes): print(esp,'Dupes:', df[date_dupes][cell_col],'\n')
        df = df[[not x for x in date_dupes]]

        df['ESP '+esp] = df[cell_col].round()
        web_df = web_df.join( df['ESP '+esp], how='outer')
    
    return web_df
    
    
def get_local_df(path):
    return pd.read_csv(path, index_col='Date',parse_dates=True)
    
    
def do_merge_df(dfA,dfB):
    return dfA.reset_index().merge(dfB.reset_index(), how='outer').set_index(dfA.index.name)


def do_default():
    this_year = datetime.now().year
    csv_file = 'data/{}.csv'.format(this_year)
    web_df = get_web_df(this_year, DEFAULT_ESPs)
    if os.path.isfile(csv_file):
        local_df = get_local_df(csv_file)
        merged_df = do_merge_df(local_df,web_df)
    else:
        merged_df = web_df
    merged_df.to_csv(csv_file)

  
if __name__ == "__main__":
    this_year = datetime.now().year
    parser = argparse.ArgumentParser(description='Pull and merge new cell counts from http://science.whoi.edu/fieldcelldata/secure/{year}/{esp}/counts.csv')
    parser.add_argument('--year', default=this_year,type=int,help='default is the current year')
    parser.add_argument('--csv-file', default='data/{}.csv'.format(this_year), help='csv file to merge new data into. default is "data/{year}.csv"')
    parser.add_argument('--ESPs', nargs='+', default=DEFAULT_ESPs, help='ESPs to pull from online. default is "{}"'.format(' '.join(DEFAULT_ESPs)))
    parser.add_argument('--output','-o', default=None, help='file to output to. default is same as --csv-file (overwrite). may also be "stdout"')
    args = parser.parse_args()
    
    web_df = get_web_df(args.year, args.ESPs)
    
    if os.path.isfile(args.csv_file):
        local_df = get_local_df(args.csv_file)
        merged_df = do_merge_df(local_df,web_df)
    else:
        merged_df = web_df
    
    if args.output == 'stdout':
        fio = StringIO()
        merged_df.to_csv(fio)
        print(fio.getvalue())
    elif args.output is None:
        merged_df.to_csv(args.csv_file)
    else:
        merged_df.to_csv(args.output)
        
        
        
        
