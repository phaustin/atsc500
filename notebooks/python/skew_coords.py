
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
from a500.skewT.skewlib import makeSkewDry
from a500.thermo.thermlib import convertTempToSkew
import datetime
import pytz

soundings= read_soundings('soundingdir')
print(soundings.keys())
print(soundings['sounding_dict'].keys())


# In[5]:


the_date=(2017,7,1,12)
the_sounding=soundings['sounding_dict'][the_date]
attributes=soundings['attributes']
#print(attributes)
fig,ax =plt.subplots(1,1,figsize=(8,8))
ax,skew = makeSkewDry(ax)
temp=the_sounding['temp']
press = the_sounding['pres']
tdew = the_sounding['dwpt']
temp_skew = convertTempToSkew(temp,press,skew)
tdew_skew = convertTempToSkew(tdew,press,skew)
ax.plot(temp_skew,press)
ax.plot(tdew_skew,press)
the_date=datetime.datetime(*the_date,tzinfo=pytz.utc)
central=pytz.timezone('US/Central')
the_date_central=the_date.astimezone(central)
title=f'Dodge City KS sounding: {str(the_date_central)}'
ax.set_title(title);
#help(convertTempToSkew)


# In[4]:


import time
print(dir(time))

