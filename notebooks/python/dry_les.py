
# coding: utf-8

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc" style="margin-top: 1em;"><ul class="toc-item"><li><span><a href="#Intro-to-netcdf" data-toc-modified-id="Intro-to-netcdf-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Intro to netcdf</a></span></li><li><span><a href="#Intro-to-python-packages" data-toc-modified-id="Intro-to-python-packages-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>Intro to python packages</a></span></li><li><span><a href="#Dumping-the-netcdf-metadata" data-toc-modified-id="Dumping-the-netcdf-metadata-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Dumping the netcdf metadata</a></span></li><li><span><a href="#problem-for-Wednesday" data-toc-modified-id="problem-for-Wednesday-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>problem for Wednesday</a></span></li></ul></div>

# In[7]:


from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np


# # Intro to netcdf
# 
# See the function descriptions and tutorial at http://unidata.github.io/netcdf4-python/

# # Intro to python packages
# 
# a. do the following to install the course python code using [pip][1]:
#    
#      cd atsc500
#      git fetch origin
#      git reset --hard origin/master
#      pip install -e .
#     
#    (this is called an "editable install", for reasons I'll explain in class)
#    
#    [1]: https://en.wikipedia.org/wiki/Pip_(package_manager)
#  
# b. Check the install by executing the cell below:
# 
#    If it succeeds, you should see:
#    
#        download case_60_10.nc: size is    499.3 Mbytes

# In[5]:


from  a500.utils.data_read import download
the_root="http://clouds.eos.ubc.ca/~phil/docs/atsc500/data/"
the_file='case_60_10.nc'
out=download(the_file,root=the_root)


# # Dumping the netcdf metadata

# Netcdf file layout:  10 groups corresponding to 10 different ensemble members.  Small slice of larger domain of LES run with surface heat flux of 60 W/m^2 and stable layer with dT/dz = 10 K/km.  Snapshots every 10 minutes for 8 hours.
# 
# We can read the metdata using 

# In[9]:


get_ipython().system('pyncdump case_60_10.nc')


# Plot $\theta$ profile for every third timestep (i.e. every 30 minutes)

# In[11]:


get_ipython().run_line_magic('matplotlib', 'inline')

def make_theta(temp,press):
    """
      temp in K
      press in Pa
      returns theta in K
    """
    p0=1.e5
    Rd=287.  #J/kg/K
    cpd=1004.  #J/kg/K
    theta=temp*(p0/press)**(Rd/cpd)
    return theta

case_name='case_60_10.nc'
#
#  look at the first ensemble member
#
ensemble='c1'
with Dataset(case_name,'r') as ncin:
    #
    # grab the group variables
    #
    group = ncin.groups['c1']
    temp=group.variables['TABS'][...]
    press=ncin.variables['press'][...]
    z=ncin.variables['z'][...]
temp=temp.mean(axis=3).mean(axis=2)

fig,ax=plt.subplots(1,1,figsize=(10,8))
for i in np.arange(0,temp.shape[0],3):
    theta = make_theta(temp[i,:],press)
    ax.plot(theta,z)
out=ax.set(xlabel=r'$\overline{\theta}$ (K)',ylabel='height (m)',
       title='LES dry run for realization 1:  surface flux=60 $W\,m^{-2}$, $\Gamma$=10 K/km')



# # problem for Wednesday
# 
# Hand in a notebook that
# 
# 1. plots the ensemble average theta profile for 1 column (1-d z vs Temp), and the ensemble averaged vertical heat flux in Watts/m^2 at z=400 meters (2-d contour plot using pcolormesh)
# 
