
# coding: utf-8

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc" style="margin-top: 1em;"><ul class="toc-item"><li><span><a href="#plot-a-sounding" data-toc-modified-id="plot-a-sounding-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>plot a sounding</a></span></li></ul></div>

# In[1]:


import numpy as np
from matplotlib import pyplot as plt
#!conda install -y pytz


# # plot a sounding

# In[2]:


from a500.soundings.wyominglib import read_soundings
from a500.skewT.fullskew import makeSkewWet,find_corners,make_default_labels
from a500.thermo.thermlib import convertTempToSkew
import datetime
import pytz

soundings= read_soundings('soundingdir')
print(soundings.keys())
print(soundings['sounding_dict'].keys())


# In[3]:


the_sounding=soundings['sounding_dict'][(2017,7,14,0)]


# In[4]:




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

skew=35.
temp=the_sounding['temp']
press = the_sounding['pres']
tdew = the_sounding['dwpt']
temp_skew = convertTempToSkew(temp,press,skew)
tdew_skew = convertTempToSkew(tdew,press,skew)

fig,ax =plt.subplots(1,1,figsize=(10,8))
corners = [10, 35]
ax, skew = makeSkewWet(ax, corners=corners, skew=skew,label_fun=label_fun)
ax.set_title('Dodge City Kansas sounding')
xcorners=find_corners(corners,skew=skew)
ax.set(xlim=xcorners,ylim=[1000,400])
ax.plot(temp_skew,press,linewidth=5)
ax.plot(tdew_skew,press,linewidth=5);

