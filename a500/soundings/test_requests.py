"""
script version of the newsoundings.ipynb notebook
"""

import requests
from bs4 import BeautifulSoup
import re
import datetime
from  datetime import timezone as tz
import pandas as pd
from collections import OrderedDict
from a405thermo.thermlib import find_esat
from a405thermo.constants  import constants as con
import numpy as np
import sys
import h5py


re_text="""
           .*Station\snumber\:\s(.+?)\n
           \s+Observation\stime\:\s(.+?)\n
           \s+Station\slatitude\:\s(.+?)\n
           \s+Station\slongitude\:\s(.+?)\n
           \s+Station\selevation\:\s(.+?)\n
           .*
        """

header_re=re.compile(re_text,re.DOTALL|re.VERBOSE)

def parse_header(header_text):
    header_text=header_text.strip()
    the_match = header_re.match(header_text)
    try:
        the_id,string_time,lat,lon,elev=the_match.groups()
    except AttributeError:
        print('parse failure with: \n',header_text)
    elev=elev.split('\n')[0]  #some soundings follow elev with Shoalwater, not Lifted
    lat=float(lat)
    lon=float(lon)
    elev=float(elev)
    day,hour = string_time.strip().split('/')
    year=int(day[:2]) + 2000
    month=int(day[2:4])
    day=int(day[4:6])
    minute=int(hour[2:])
    hour=int(hour[:2])
    theDate=datetime.datetime(year,month,day,hour,minute,0,0,tz.utc)
    return theDate,the_id,lat,lon,elev

def parse_data(data_text):
    all_lines=data_text.strip().split('\n')
    count=0
    theLine=all_lines[count]
    try:
        while theLine.find('PRES   HGHT   TEMP   DWPT') < 0:
            count += 1
            theLine = all_lines[count]
        header_names=all_lines[count].lower().split()
    except IndexError:
        print("no column header line found in sounding")
        sys.exit(1)
    count += 1  #go to unit names
    unit_names=all_lines[count].split()
    count+=2  #skip a row of ------
    data_list=[]
    while True:
        try:
            the_line=all_lines[count]
            dataFields = the_line.split()
            if len(dataFields) == 11:
                try:
                    dataFields = [float(number) for number in dataFields]
                    es = find_esat(dataFields[3] + 273.15)*0.01  #get vapor pressure from dewpoint in hPa
                    dataFields[5] = (con.eps*es/(dataFields[0] - es))*1.e3   #g/kg
                except VauleError:
                    print('trouble converting dataFields to float')
                    print(dataFields)
                    sys.exit(1)
                data_list.append(dataFields)
                #
                # get the next line
                #
            count += 1
            theLine = all_lines[count]
        except IndexError:
            break
    df_out=pd.DataFrame.from_records(data_list,columns=header_names)
    return df_out,unit_names

def make_frames(html_doc):
    """
      input: web page from wyoming upperair sounding site
             http://weather.uwyo.edu/cgi-bin/sounding retrieved by
             the requests module
      output: attr_dict dictionary with ['header', 'site_id','longitude','latitude', 'elevation', 'units']
              sound_dict soudning dictionary with sounding times as keys and sounding as dataframes
    """  
    soup=BeautifulSoup(html_doc,'html.parser')
    keep=list()
    header= (soup.find_all('h2'))[0].text
    print('header is: ',header)
    for item in soup.find_all('pre'):
        keep.append(item.text)

    sounding_dict=OrderedDict()
    evens=np.arange(0,len(keep),2)
    for count in evens:
        df,units=parse_data(keep[count])
        theDate,the_id,lat,lon,elev=parse_header(keep[count+1])    
        sounding_dict[theDate]=df
        
    attr_dict=dict(units=';'.join(units),site_id=the_id,
                   latitude=lat,longitude=lon,elevation=elev,
                   header = header)
    
    return attr_dict,sounding_dict
    


if __name__ == "__main__":


    url_template=("http://weather.uwyo.edu/cgi-bin/sounding?"
              "region={region:s}"
              "&TYPE=TEXT%3ALIST"
              "&YEAR={year:s}"
              "&MONTH={month:s}"
              "&FROM={start:s}"
              "&TO={stop:s}"
              "&STNM={station:s}")

    values=dict(region='samer',year='2013',month='2',start='0100',stop='2800',station='82965')
    #values=dict(region='nz',year='2013',month='2',start='0100',stop='2800',station='93417')
    #values=dict(region='naconf',year='2013',month='2',start='0100',stop='2800',station='71802')
    #values=dict(region='ant',year='2013',month='07',start='0100',stop='2800',station='89009')
    values=dict(region='naconf',year='2012',month='3',start='0100',stop='0212',station='72340')
    #naconf, samer, pac, nz, ant, np, europe,africa, seasia, mideast
    url=url_template.format_map(values)

    do_web = True
    backup_file='backup.txt'
    if do_web:
        html_doc = requests.get(url).text
        print(len(html_doc))
        with open(backup_file,'w') as f:
            f.write(html_doc)
        if len(html_doc) < 2000:
            print('debug: short html_doc, something went wrong:',html_doc)
            sys.exit(1)
    else:
        with open(backup_file,'r') as f:
            html_doc=f.read()

    attr_dict,sounding_dict = make_frames(html_doc)
    attr_dict['history']="written by test_requests.py"
    key_list=['header', 'site_id','longitude','latitude', 'elevation', 'units','history']

    name = 'single_littlerock.h5'    
    with pd.HDFStore(name,'w') as store:
        for key,value in sounding_dict.items():
            thetime=key.strftime("Y%Y_%b_%d_%HZ")
            store.put(thetime,value,format='table')

    with h5py.File(name,'a') as f:
        for key in key_list:
            print('writing key, value: ',key,attr_dict[key])
            f.attrs[key]=attr_dict[key]
        f.close()

    with h5py.File(name,'r') as f:
        keys=f.attrs.keys()
        for key in keys:
            try:
                print(key,f.attrs[key])
            except OSError:
                pass





