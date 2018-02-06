
# coding: utf-8

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc" style="margin-top: 1em;"><ul class="toc-item"><li><span><a href="#Supress-autoscrolling" data-toc-modified-id="Supress-autoscrolling-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Supress autoscrolling</a></span><ul class="toc-item"><li><span><a href="#Grab-a-Little-Rock-sounding" data-toc-modified-id="Grab-a-Little-Rock-sounding-1.1"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>Grab a Little Rock sounding</a></span></li><li><span><a href="#Select-one-sounding" data-toc-modified-id="Select-one-sounding-1.2"><span class="toc-item-num">1.2&nbsp;&nbsp;</span>Select one sounding</a></span></li><li><span><a href="#Save-the-metadata-for-plotting" data-toc-modified-id="Save-the-metadata-for-plotting-1.3"><span class="toc-item-num">1.3&nbsp;&nbsp;</span>Save the metadata for plotting</a></span></li><li><span><a href="#Convert-temperature-and-dewpoint-to-skew-coords" data-toc-modified-id="Convert-temperature-and-dewpoint-to-skew-coords-1.4"><span class="toc-item-num">1.4&nbsp;&nbsp;</span>Convert temperature and dewpoint to skew coords</a></span></li><li><span><a href="#Plot-the-sounding,-making-the-sounding-lines-thicker" data-toc-modified-id="Plot-the-sounding,-making-the-sounding-lines-thicker-1.5"><span class="toc-item-num">1.5&nbsp;&nbsp;</span>Plot the sounding, making the sounding lines thicker</a></span></li><li><span><a href="#turn-off-log(0)-warning" data-toc-modified-id="turn-off-log(0)-warning-1.6"><span class="toc-item-num">1.6&nbsp;&nbsp;</span>turn off log(0) warning</a></span></li><li><span><a href="#find-the-$\theta_{es}$-of-the--LCL" data-toc-modified-id="find-the-$\theta_{es}$-of-the--LCL-1.7"><span class="toc-item-num">1.7&nbsp;&nbsp;</span>find the $\theta_{es}$ of the  LCL</a></span><ul class="toc-item"><li><span><a href="#What-is-the-LCL-of-this-air?" data-toc-modified-id="What-is-the-LCL-of-this-air?-1.7.1"><span class="toc-item-num">1.7.1&nbsp;&nbsp;</span>What is the LCL of this air?</a></span></li></ul></li><li><span><a href="#Find-the-vertical-buoyancy-profile" data-toc-modified-id="Find-the-vertical-buoyancy-profile-1.8"><span class="toc-item-num">1.8&nbsp;&nbsp;</span>Find the vertical buoyancy profile</a></span></li><li><span><a href="#get-the-limits-of-integration" data-toc-modified-id="get-the-limits-of-integration-1.9"><span class="toc-item-num">1.9&nbsp;&nbsp;</span>get the limits of integration</a></span></li><li><span><a href="#Calculate-the-cumulative-cape" data-toc-modified-id="Calculate-the-cumulative-cape-1.10"><span class="toc-item-num">1.10&nbsp;&nbsp;</span>Calculate the cumulative cape</a></span></li><li><span><a href="#interpolating-between-z-levels" data-toc-modified-id="interpolating-between-z-levels-1.11"><span class="toc-item-num">1.11&nbsp;&nbsp;</span>interpolating between z levels</a></span></li></ul></li></ul></div>

# # Supress autoscrolling

# In[1]:


get_ipython().run_cell_magic('javascript', '', 'IPython.OutputArea.prototype._should_scroll = function(lines) {\n    return false;\n}')


# In[2]:


import numpy as np
import pandas as pd
from pprint import pformat

from a405.thermo.constants import constants as c
from a405.thermo.thermlib import convertSkewToTemp, convertTempToSkew
from a405.skewT.fullskew import makeSkewWet,find_corners,make_default_labels


# In[3]:


from a405.soundings.wyominglib import write_soundings, read_soundings
from matplotlib import pyplot as plt


# ## Grab a Little Rock sounding
# 
# Set get_data = True to fetch the Wyoming sounding data from the website and store it.

# In[4]:


get_data = False
values=dict(region='naconf',year='2012',month='7',start='0100',stop='3000',station='72340')
if get_data:
    write_soundings(values, 'littlerock')
soundings= read_soundings('littlerock')


# In[5]:


soundings['sounding_dict'].keys()


# ## Select one sounding

# In[6]:


the_time=(2012,7,17,0)
sounding=soundings['sounding_dict'][the_time]
sounding.columns


# ## Save the metadata for plotting

# In[7]:


title_string=soundings['attributes']['header']
index=title_string.find(' Observations at')
location=title_string[:index]
print(f'location: {location}')

units=soundings['attributes']['units'].split(';')
units_dict={}
for count,var in enumerate(sounding.columns[1:]):
    units_dict[var]=units[count]
#
# use the pretty printer to print the dictionary
#
print(f'units: {pformat(units_dict)}')


# ## Convert temperature and dewpoint to skew coords

# In[8]:


skew=35.
triplets=zip(sounding['temp'],sounding['dwpt'],sounding['pres'])
xcoord_T=[]
xcoord_Td=[]
for a_temp,a_dew,a_pres in triplets:
    xcoord_T.append(convertTempToSkew(a_temp,a_pres,skew))
    xcoord_Td.append(convertTempToSkew(a_dew,a_pres,skew))


# ## Plot the sounding, making the sounding lines thicker

# In[9]:


def label_fun():
    """
    override the default rs labels with a tighter mesh
    """
    from numpy import arange
    #
    # get the default labels
    #
    tempLabels,rsLabels, thetaLabels, thetaeLabels = make_default_labels()
    #
    # change the temperature and rs grids
    #
    tempLabels = range(-40, 50, 2)
    rsLabels = [0.1, 0.25, 0.5, 1, 2, 3] + list(np.arange(4, 28, 2)) 
    return tempLabels,rsLabels, thetaLabels, thetaeLabels


# In[10]:


fig,ax =plt.subplots(1,1,figsize=(8,8))
corners = [10, 35]
ax, skew = makeSkewWet(ax, corners=corners, skew=skew,label_fun=label_fun)
#ax,skew = makeSkewWet(ax,corners=corners,skew=skew)
out=ax.set(title=title_string)
xcorners=find_corners(corners,skew=skew)
ax.set(xlim=xcorners,ylim=[1000,400]);
l1,=ax.plot(xcoord_T,sounding['pres'],color='k',label='temp')
l2,=ax.plot(xcoord_Td,sounding['pres'],color='g',label='dew')
[line.set(linewidth=3) for line in [l1,l2]];


# ## turn off log(0) warning

# In[11]:


np.seterr(all='ignore');


# ## find the $\theta_{es}$ of the  LCL

# In[12]:


from a405.thermo.thermlib import find_Tmoist,find_thetaep,find_rsat,find_Tv,find_lcl,find_thetaes,find_thetaet
#
# find thetae of the surface air, at index 0
#
sfc_press,sfc_temp,sfc_td =[sounding[key][0] for key in ['pres','temp','dwpt']]
#
sfc_press,sfc_temp,sfc_td = sfc_press*100.,sfc_temp+c.Tc,sfc_td+c.Tc


# ### What is the LCL of this air?

# In[13]:


Tlcl, plcl=find_lcl(sfc_td, sfc_temp,sfc_press)


# In[14]:


print(f'found Tlcl={Tlcl} K, plcl={plcl} Pa')


# In[15]:


#  convert to mks and find surface rv and thetae
#
sfc_rvap = find_rsat(sfc_td,sfc_press)
lcl_rvap = find_rsat(Tlcl,plcl)
sfc_thetae=find_thetaes(Tlcl,plcl)
press=sounding['pres'].values*100.
#
# find the index for 70 hPa pressure -- searchsorted requires
# the pressure array to be increasing, so flip it for the search,
# then flip the index.  Above 70 hPa thetae goes bananas, so
# so trim so we only have good values
#
toplim=len(press) - np.searchsorted(press[::-1],.7e4)
clipped_press=press[:toplim]
#
# find temps, rv and rl along that adiabat
#
adia_temps= np.array([find_Tmoist(sfc_thetae,the_press) 
                      for the_press in clipped_press])
adia_rvaps = find_rsat(adia_temps,clipped_press)
adia_rls = sfc_rvap - adia_rvaps
env_temps = (sounding['temp'].values + c.Tc)[:toplim]
env_Td = (sounding['dwpt'].values + c.Tc)[:toplim]
height = sounding['hght'].values[:toplim]
pairs = zip(env_Td,clipped_press)
env_rvaps= np.array([find_rsat(td,the_press) for td,the_press in pairs])
env_Tv = find_Tv(env_temps,env_rvaps)
adia_Tv = find_Tv(adia_temps,adia_rvaps,adia_rls)
xcoord_thetae=[]
press_hPa = clipped_press*1.e-2
#
# convert the adiabatic thetae sounding to skewT coords
#
for a_temp,a_press in zip(adia_temps - c.Tc,press_hPa):
    out=convertTempToSkew(a_temp,a_press,skew)
    xcoord_thetae.append(out)
ax.plot(xcoord_thetae,press_hPa,color='r',label='thetae',linewidth=3.)
xlcl=convertTempToSkew(Tlcl - c.Tc,plcl*0.01,skew)
ax.plot(xlcl,plcl*0.01,'bo',markersize=8)
display(fig)


# ## Find the vertical buoyancy profile
# 
# The cloud parcel will begin decelerating when the buoyancy becomes negative.  It will overshoot that level and rise into the inversion.

# In[16]:


def find_buoy(adia_Tv,env_Tv):
  buoy=c.g0*(adia_Tv - env_Tv)/env_Tv
  return buoy

#
# moved find_buoy to library
#
from a405.thermo.thermlib import find_buoy
    
fig,ax =plt.subplots(1,1)
buoy=find_buoy(adia_Tv,env_Tv)
ax.plot(buoy,clipped_press*1.e-2)
ax.invert_yaxis()
out=ax.set(ylabel='press (hPa)',xlabel='buoyancy (m/s^2)')
ax.grid(which='both')


# ## get the limits of integration

# In[17]:


#np.searchsorted finds the first crossing
zerobot=np.searchsorted(buoy,0)
#flip the array and search backwards for top crossing
zerotop=len(buoy) - np.searchsorted(buoy[::-1],0)
print(zerotop, toplim)
print('pressure levels for crossing: {} hPa {} hPa'      .format(press[zerobot]*0.01,press[zerotop]*0.01))


# In[18]:


clipped_buoy = buoy[zerobot:zerotop]
clipped_height = height[zerobot:zerotop]
#average the levels to get layer buoyancy
layer_buoy = (clipped_buoy[:-1] + clipped_buoy[1:])/2.
cape = np.sum(layer_buoy*np.diff(clipped_height))
print('cape is {:6.2f} J/kg'.format(cape))


# ## Calculate the cumulative cape
# 
# We want to see how CAPE accumulates as the parcel ascends to find the level in the inversion at which CAPE returns to zero, i.e. the maximum cloud top. Do this integral in python
# 
# $$CAPE(z) = \int_{zbot}^{z} B dz$$

# In[19]:


clipped_buoy = buoy[zerobot:]
clipped_height = height[zerobot:]
new_press= clipped_press[zerobot:]
#average the levels to get layer buoyancy
layer_buoy = (clipped_buoy[:-1] + clipped_buoy[1:])/2.
layer_height = (clipped_height[:-1] + clipped_height[:-1])/2.
layer_press = (new_press[:-1] + new_press[1:])/2.
cape = np.cumsum(layer_buoy*np.diff(clipped_height))
fig,ax=plt.subplots(1,1,figsize=(8,8))
ax.plot(cape,layer_press*1.e-2)
ax.grid(True,which='both')
ax.set(xlabel='cumulative cape ($m^2/s^2$)',ylabel='pressure (hPa)')
ax.invert_yaxis()


# ## interpolating between z levels
# 
# We want to be able to solve for the buoyancy at any height so we can use a ordinary differential equations to describe
# mixing.   The scipy module has a a one-dimensional interpolator interp1d that we can use to interpolate between
# the pressure levels reported by the balloon

# In[20]:


from scipy.interpolate import interp1d
env_Tv_interp=interp1d(clipped_press,env_Tv)
fig,ax = plt.subplots(1,1)
ax.plot(env_Tv,clipped_press*0.01,label='orig sounding')
ax.invert_yaxis()
reg_press=np.linspace(9.8e4,6.e4,25)
ax.plot(env_Tv_interp(reg_press),reg_press*0.01,'ro',label='interp')
ax.set(ylim=(1000,650),xlim=(260,320))
ax.grid(True)
out=ax.legend(loc='best')

