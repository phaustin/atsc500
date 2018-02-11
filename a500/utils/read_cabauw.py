import glob
from netCDF4 import Dataset
from dateutil.parser import parse
import datetime
import numpy as np
import matplotlib.dates as mdates
from a500.utils.data_read import download
from pathlib import Path
from collections import defaultdict
import pdb
from netCDF4 import Dataset
import os


def date_sort(filepath):
    """
    given a filepath, return a tuple with the year and month for a cabauw nc file 
    i.e. given cesar_tower_meteo_lc1_t10_v1.0_201401.nc return (2014,1)

    Parameters
    ----------

    filepath: str or Path object

    Returns
    -------

    year, month: tuple of ints
      e.g.  (2014,10)
    
    """
    filestr=str(filepath)
    #split string like: _201401.nc
    int_date=(filestr.split('_')[-1]).split('.nc')[0]
    year = int(int_date[:4])
    month = int(int_date[4:6])
    return year,month


def make_date(ncfile):
    """
    open a cabauw file and turn the time information which is in the form
    product:date_start_of_data = "2014-01-01T00:00Z" ; plus
    time:units = "hours since 2014-01-01 00:00:00 0:00" ;

    into a numpy array of python datetime objects

    Parameters
    ----------

    ncfile:  open netcdf Dataset object

    Returns
    -------

    time_vec: ndarray
       array of datetime objects with shape [num_days,24,6]
       where num_days is the nunmber of days in the month
       24 is the number of hours in a day and 6 is the number
       of 10 minute measurements in an hour
    """
    the_time=ncfile.variables['time'][...]
    start_date=ncfile.variables['product'].date_start_of_data
    start_date = parse(start_date)
    time_vec=[]
    for the_hour in the_time:
        time_vec.append(start_date + datetime.timedelta(hours=float(the_hour)))
    time_vec=np.array(time_vec)
    time_vec=time_vec.reshape(-1,24,6)
    return time_vec


def make_file_dict(root_dir,year):
    """
    return a dictionary of file paths keyed by cabauw file type, which can be
       'surface_fluxes', 'surface meteorological' or 'tower_meteorological'

    Parameters
    ----------

    root_dir: str
      path to folder holding the cabauw netcdf files

    year: int
      4 digit integer year, e.g. 2014


    Returns
    -------

    file_dict: dict
       dictionary with 3 keys, each key returns a list of filenames as Path objects
    

    """
    file_dict = {'surface_fluxes':f"**/*surface_flux*{year}??*.nc",
                'surface_meteorological':f"**/*surface_meteo*{year}??*.nc",
                'tower_meteorological':f"**/*tower_meteo*{year}??*.nc"}
    for file_type,file_str in file_dict.items():
        file_list=list(root_dir.glob(file_str))
        file_list.sort(key=date_sort)
        file_dict[file_type]=file_list
    return file_dict

def get_attrs(ncvar):
    """

    Utilicty function to return all attributes for a netcdf variable as
    a dictionray keyed by the variable name

    Parameters
    ----------

    ncvar: open netcdf4 variable object

    Returns
    -------

    attr_dict: dict
       dictionary with key,item pairs for each attribute
    """
    attributes=ncvar.ncattrs()
    attr_dict={}
    for attr in attributes:
        attr_dict[attr]=getattr(ncvar,attr)
    return attr_dict

def write_dates(filetype,filelist):
    """
    give a list of files of filetype 
    'surface_fluxes', 'surface meteorological' or 'tower_meteorological'
    loop over the filepaths in filelist and extract a dictionary with basic
    time and location data for that filelist

    Parameters
    ----------

    filetype: str
       one of 'surface_fluxes', 'surface meteorological' or 'tower_meteorological'
    filelist: list of str or Path objects with path to file

    Returns
    -------

    data_dict: dict
        nested dictionary keyed by [filetype](year,month) e.g. (2014,10) with each entry
        another dictionary with keys 

        name=str_file,timevec=the_time,lon=lon,lat=lat,
            filetype=filetype,start_date=start_date,attributes=attr_dict
    
    """
    data_dict=defaultdict(dict)
    for the_file in filelist:
        str_file=str(the_file)
        with Dataset(str_file,'r') as ncin:
            details=ncin.variables['iso_dataset']
            attr_dict=get_attrs(details)
            lon=attr_dict['westbound_longitude']
            lat=attr_dict['northbound_latitude']
            the_time = make_date(ncin)
            start_date=ncin.variables['product'].date_start_of_data
            start_date=parse(start_date)
            start_month=date_sort(str_file)
        data_dict[start_month]=dict(name=str_file,timevec=the_time,lon=lon,lat=lat,
                                    filetype=filetype,start_date=start_date)
    return data_dict

def store_months(data_dict):
    """
    given a nested dictionary with keys data_dict[filetype][year,month]
    where filetype is one of 
    'surface_fluxes', 'surface meteorological' or 'tower_meteorological'
    and [year,month] is a key like (2014,10) points to a dictionary returned
    by write_dates  get fluxes and tower data for the monthly files from cabauw
    """
    full_dict=defaultdict(dict)
    for filetype in data_dict.keys():
        var_dict=data_dict[filetype]
        month_list=list(var_dict.keys())
        out_dict=defaultdict(dict)
        var_attrs=defaultdict(dict)
        full_dict[filetype]['filelist']=[]
        for the_month in month_list:
            month_dict=var_dict[the_month]
            full_dict[filetype]['filelist'].append(month_dict['name'])
            print('working on {filetype:}: {start_date:} '.format_map(month_dict))
            if filetype == 'tower_meteorological':
                with Dataset(month_dict['name'],'r') as ncin:
                    if 'has_time' not in out_dict[the_month]:
                        out_dict[the_month]['timevec'] = month_dict['timevec']
                        out_dict[the_month]['has_time']=True
                    out_dict[the_month]['z'] = ncin.variables['z'][...]
                    for var in ['F','TA','TD','Q','D']:
                        out_dict[the_month][var] = ncin.variables[var][...].reshape(-1,24,6,7)
                        var_attrs[var]=get_attrs(ncin.variables[var])
            elif filetype == 'surface_fluxes': 
                with Dataset(month_dict['name'],'r') as ncin:
                    for var in ['H','UST','LE']:
                        out_dict[the_month][var] = ncin.variables[var][...].reshape(-1,24,6)
                        var_attrs[var]=get_attrs(ncin.variables[var])
            elif filetype == 'surface_meteorological': 
                with Dataset(month_dict['name'],'r') as ncin:
                    for var in ['P0','TA002','Q002','F010']:
                        out_dict[the_month][var] = ncin.variables[var][...].reshape(-1,24,6)
                        var_attrs[var]=get_attrs(ncin.variables[var])
            else:
                raise ValueError(f"didn't recognize {filetype}")
        full_dict[filetype]['data']=out_dict
        full_dict[filetype]['var_attrs']=var_attrs
        full_dict[filetype]['lat']=var_dict[the_month]['lat']
        full_dict[filetype]['lon']=var_dict[the_month]['lon']
    return full_dict


# We'll need to know sunrise, sunset and solar noon to interpret our
# data.  Here is how you find these with the 
# [pyephem](http://stackoverflow.com/questions/2637293/calculating-dawn-and-sunset-times-using-pyephem) module

# In[8]:

def sunset_sunrise():
    import ephem
    for the_month in month_dict.keys():
        var='tower_meteorological'
        start_time=month_dict[the_month]['timevec'][0,0,0]
        cabauw=ephem.Observer()
        cabauw.date=start_time
        cabauw.lon = var_attrs['lon']
        cabauw.lat = var_attrs['lat']
        sunrise=cabauw.next_rising(ephem.Sun())
        noon = cabauw.next_transit(ephem.Sun(),start=sunrise)
        sunset = cabauw.next_setting(ephem.Sun())
        print('sunrise is {} UTC'.format(sunrise))
        print('solar noon {} UTC'.format(noon))
        print('sunset is {} UTC'.format(sunset))


# # Data first look

# In[9]:


# from matplotlib import pyplot as plt
# plt.style.use('ggplot')
# plt.close('all')
# print('starting: ',month_dict.keys())    
# for the_month in month_dict.keys():
#     hourly_wind_avg = month_dict[the_month]['F'].mean(axis=2)
#     z=month_dict[the_month]['z']
#     hour=2
#     fig,ax=plt.subplots(1,2,figsize=(8,6))
#     ax[0].plot(hourly_wind_avg[:,hour,:].T,z)
#     ax[0].set(title='hour: {} UTC'.format(hour))
#     hour=14
#     ax[1].plot(hourly_wind_avg[:,hour,:].T,z)
#     ax[1].set(title='hour: {} UTC'.format(hour))
#     fig.suptitle('{} hourly avg winds'.format(the_month))

#     #
#     # date plotting tips at http://matplotlib.org/users/recipes.html
#     #
#     the_time = month_dict[the_month]['timevec']
#     H=month_dict[the_month]['H']
#     fig,ax=plt.subplots(1,1,figsize=(8,6))
#     fig.autofmt_xdate()
#     ax.plot(the_time.flat,H.flat)
#     title='sensible heat flux for {}'.format(the_month)
#     ax.set(title=title,ylabel='H $(W\,m^{-2})$')
# print('finished plot: ',month_dict.keys())


# # Checkpoint all the data into one netcdf file
# 
# I want to save month_dict into a netcdf file so we don't need
# to repeat this processing but can start with a merged dataset
# that has all days of interest and all
# instruments in a single place.   To do that, I group the measurements
# into individual days using [netcdf groups](http://unidata.github.io/netcdf4-python)
# 
# I transfer all the attributes I read into the var_attrs dict so I maintain
# the original metadata as much as possible

# In[10]:


def write_attrs(ncvar,attr_dict):
    for name,item in attr_dict.items():
        if name != '_FillValue':
            setattr(ncvar,name,item)
    return None

def write_var(full_dict,filetype,date_group,the_month,days_name,varname,dim_tup):
        the_data=full_dict[filetype]['data'][the_month][varname]
        var_nc=date_group.createVariable(varname,the_data.dtype,dim_tup)
        var_nc[...]=the_data[...]
        attrs=full_dict[filetype]['var_attrs'][varname]
        write_attrs(var_nc,attrs)


def write_cdf(full_dict,outfile='caubau_ubc.nc'):
    if os.path.exists(outfile):
        os.remove(outfile)   
    with Dataset(outfile,'w') as ncout:
        the_months=list(full_dict['surface_fluxes']['data'].keys())
        for the_month in the_months:
            speed=full_dict['tower_meteorological']['data'][the_month]['F']
            days_name='days{:2d}'.format(speed.shape[0])
            dimnames=[days_name,'hours','min10','z']
            dim_info=zip(dimnames,speed.shape)
            for name,length in dim_info:     
                try:
                    ncout.createDimension(name,length) 
                except RuntimeError:
                    #
                    # ignore if dimensions already exist
                    #
                    pass
            nc_month='m{0:02d}_{1:02d}'.format(*the_month)
            date_group=ncout.createGroup(nc_month)
            setattr(date_group,'month',np.array(the_month,dtype=np.int32))
            dim_tup=[days_name,'hours','min10']
            for var in ['H','LE','UST']:
                filetype='surface_fluxes'
                write_var(full_dict,filetype,date_group,the_month,days_name,var,dim_tup)
            for var in ['P0','TA002','Q002','F010']:
                filetype='surface_meteorological'
                write_var(full_dict,filetype,date_group,the_month,days_name,var,dim_tup)
            dim_tup=[days_name,'hours','min10','z']
            for var in ['TA','D','Q','TD','F']:
                filetype='tower_meteorological'
                write_var(full_dict,filetype,date_group,the_month,days_name,var,dim_tup)
            the_time=full_dict[filetype]['data'][the_month]['timevec']
            float_time=np.array([item.timestamp() for item in the_time.flat])
            float_time=float_time.reshape(-1,24,6)   
            time_nc=date_group.createVariable('time',float_time.dtype,
                                         [days_name,'hours','min10']) 
            time_nc[...]=float_time[...]
            time_nc.timezone='UTC'
            time_nc.comment='convert using datetime.fromtimestamp(var,pytz.utc)'
            try:
                filetype = 'tower_meteorological'
                z=full_dict[filetype]['data'][the_month]['z']
                z_nc=ncout.createVariable('z',z.dtype,['z'])
                z_nc[...]=z[...]
            except RuntimeError:
                #
                # only create ncout z vector once
                #
                pass
            ncout.history='written by read_cabauw.ipynb'
            #full_dict[filetype]['var_attrs']
            for var in ['lat','lon']:
                setattr(ncout,var,float(full_dict[filetype][var]))
            ncout.lat_units='degrees north'
            ncout.lon_units='degrees east'
            for filetype in full_dict.keys():
                attrname=f'{filetype}_filelist'
                def split_file(filename):
                    filename=Path(filename)
                    filename=str(filename.parts[-1])
                    return filename
                filelist = [split_file(item) for item in full_dict[filetype]['filelist']]
                filelist=','.join(filelist)
                setattr(ncout,attrname,filelist)

def main():
    root_dir = "/Users/phil/Downloads/data"
    root_dir=Path(root_dir)
    year = 2014
    file_dict=make_file_dict(root_dir,year)
    keep_months={}
    for filetype,file_list in file_dict.items():
        keep_months[filetype]=write_dates(filetype,file_list)
    print(list(keep_months.keys()))
    full_dict=store_months(keep_months)
    write_cdf(full_dict)


if __name__ == "__main__":
    output=main()




