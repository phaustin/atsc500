
# coding: utf-8

# In[1]:


from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np


# In[2]:


from  a500.utils.data_read import download
the_root="http://clouds.eos.ubc.ca/~phil/docs/atsc500/data/"
the_file='case_60_10.nc'
out=download(the_file,root=the_root)


# # Working with groups
# 
# the Dataset object has a method call groups that returns a dictionary with all groups
# 
# Similarly, each group object has a method called variables that returns a dictionary of all variables
# 
# So to get the shape of the wvel array for each group member, do something like this:

# In[3]:


case_name='case_60_10.nc'
#
# get the names of all the groups in the file
#
with Dataset(case_name,'r') as ncin:
    #
    # grab the group variables
    #
    groupnames=list(ncin.groups.keys())
    print(groupnames)
#
# iterate over the groups dictionary and get vvel
#
with Dataset(case_name,'r') as ncin:
    for name,group in ncin.groups.items():
        vvel = group.variables['W']
        print(vvel.shape)

