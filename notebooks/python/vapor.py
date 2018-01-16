
# coding: utf-8

# # Stull Ch. 2 problem 8

# In[10]:


cpd=1004.
rho=1
h_flux_bot=rho*0.2*cpd
h_flux_top=rho*0.1*cpd
net_heating=h_flux_bot-h_flux_top  #W/m^2
delta_z=100  #meters
total_mass=rho*delta_z  #kg/m^2
dT_dt=net_heating/(total_mass*cpd)*3600 #K/hour  
print('net heating is {:5.1f} K/hour'.format(dT_dt))


# # Getting started with vapor

# 1.  Fetch the netcdf file tropical.nc from my webserver (50 Mbytes)

# In[2]:


from a500.utils.data_read import download
the_root = "http://clouds.eos.ubc.ca/~phil/docs/atsc500"
the_file =  "tropical_subset.nc"
out=download(the_file,root=the_root)


# 2.  Identify the file structure using ncdump

# In[9]:


get_ipython().system('pyncdump tropical_subset.nc')


# 3.  Read a variable using the netCDF4 module:  http://unidata.github.io/netcdf4-python/

# In[4]:


from netCDF4 import Dataset
#
# the variables attribute returns a dictionary of variables
#
with Dataset('tropical_subset.nc','r') as nc_in:
    print(list(nc_in.variables.keys()))
    the_temp=nc_in.variables['TABS'][...]
    print(the_temp.shape)
    #
    # remove the time dimension, since we only have one timestep
    #
    the_temp=the_temp.squeeze()
    print(the_temp.shape)
    the_height=nc_in.variables['z'][...]
    print(the_height.shape)


# #Plotting three temperature profiles

# In[5]:


get_ipython().run_line_magic('matplotlib', 'inline')
#
# plot 3 arbitrary columns
#
from matplotlib import pyplot as plt
meter2km=1.e-3
random_xy=[(10,20),(80,40),(25,75)]
fig,ax=plt.subplots(1,1)
for x,y in random_xy:
    temp_profile=the_temp[:,x,y]
    ax.plot(temp_profile,the_height*meter2km)  
ax.set(xlabel='temperature (K)',ylabel='height (km)',title='3 profiles')


# # Calculate the horizontal average temperature

# In[6]:


print(the_temp.shape)
xmean=the_temp.mean(axis=2)
print(xmean.shape)
xymean=xmean.mean(axis=1)
print(xymean.shape)


# ## Calculate the temperature pertubation

# In[7]:


fig,ax=plt.subplots(1,1)
for x,y in random_xy:
    temp_profile=the_temp[:,x,y] - xymean
    ax.plot(temp_profile,the_height*meter2km)  
ax.set(xlabel='temperature (K)',ylabel='height (km)',title='3 perturbation profiles')


# # Getting started with vapor

# 1. Install vapor per [these instructions](https://www.vapor.ucar.edu/docs/usage/getting-started-vapor)
# 1. Download and unzip the following file:

# In[8]:


root = "http://clouds.eos.ubc.ca/~phil/docs/atsc500" 
file = "tropical.zip"
download(file,root=root)


# Start vapor, go to "Data -> Load dataset into current session" and see if you get a
# 3-d visualization of TABS
