
# coding: utf-8

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc" style="margin-top: 1em;"><ul class="toc-item"><li><span><a href="#Scale-heights-for-typical-atmospheric-soundings" data-toc-modified-id="Scale-heights-for-typical-atmospheric-soundings-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Scale heights for typical atmospheric soundings</a></span><ul class="toc-item"><li><span><a href="#Plot-McClatchey's-US-Standard-Atmospheres" data-toc-modified-id="Plot-McClatchey's-US-Standard-Atmospheres-1.1"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>Plot McClatchey's US Standard Atmospheres</a></span></li></ul></li><li><span><a href="#McClatchey-US-standard-atmospheres" data-toc-modified-id="McClatchey-US-standard-atmospheres-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>McClatchey US standard atmospheres</a></span><ul class="toc-item"><li><span><a href="#Data-layout" data-toc-modified-id="Data-layout-2.1"><span class="toc-item-num">2.1&nbsp;&nbsp;</span>Data layout</a></span><ul class="toc-item"><li><span><a href="#Plot--temp-and-vapor-mixing-ratio-rmix-($\rho_{H2O}/\rho_{air}$)" data-toc-modified-id="Plot--temp-and-vapor-mixing-ratio-rmix-($\rho_{H2O}/\rho_{air}$)-2.1.1"><span class="toc-item-num">2.1.1&nbsp;&nbsp;</span>Plot  temp and vapor mixing ratio rmix ($\rho_{H2O}/\rho_{air}$)</a></span></li></ul></li><li><span><a href="#Calculating-scale-heights-for-temperature-and-air-density" data-toc-modified-id="Calculating-scale-heights-for-temperature-and-air-density-2.2"><span class="toc-item-num">2.2&nbsp;&nbsp;</span>Calculating scale heights for temperature and air density</a></span><ul class="toc-item"><li><span><a href="#How-do-$\overline{H_p}$-and-$\overline{H_\rho}$-compare-for-the-tropical-sounding?" data-toc-modified-id="How-do-$\overline{H_p}$-and-$\overline{H_\rho}$-compare-for-the-tropical-sounding?-2.2.1"><span class="toc-item-num">2.2.1&nbsp;&nbsp;</span>How do $\overline{H_p}$ and $\overline{H_\rho}$ compare for the tropical sounding?</a></span></li><li><span><a href="#How-well-do-these-average-values-represent-the-pressure-and-density-profiles?" data-toc-modified-id="How-well-do-these-average-values-represent-the-pressure-and-density-profiles?-2.2.2"><span class="toc-item-num">2.2.2&nbsp;&nbsp;</span>How well do these average values represent the pressure and density profiles?</a></span></li><li><span><a href="#Assignment-for-Friday" data-toc-modified-id="Assignment-for-Friday-2.2.3"><span class="toc-item-num">2.2.3&nbsp;&nbsp;</span>Assignment for Friday</a></span></li></ul></li></ul></li></ul></div>

# # Scale heights for typical atmospheric soundings

# ## Plot McClatchey's US Standard Atmospheres
# 
# There are five different average profiles for the tropics, subtropical summer, subtropical winder, midlatitude summer, midlatitude winter.  These are called the US Standard Atmospheres.  This notebook shows how to read and plot the soundings, and calculate the pressure and density scale heights.

# In[1]:


from matplotlib import pyplot as plt
import matplotlib.ticker as ticks
import numpy as np
from a500.utils.data_read import download
from pathlib import Path
import json
import pandas as pd
from zipfile import ZipFile


# # McClatchey US standard atmospheres

# Download std_soundings.zip, which is a folder containing 5 standard atmosperic soundings
# 
# Extract the folder from the zipfile and read the csv files into pandas dataframes.  To learn
# more about dataframes, see  the [Jake Vanderplas pandas chapter](https://github.com/jakevdp/PythonDataScienceHandbook/blob/599aa0fe3f882c0001670e676e5a8d43b92c35fc/notebooks/03.01-Introducing-Pandas-Objects.ipynb)

# In[2]:


the_root="http://clouds.eos.ubc.ca/~phil/docs/atsc500/data"
the_file='std_soundings.zip'
out=download(the_file,root=the_root)
#
# use pathlib.Path to manipulate the unix/windows file paths
#
the_file = Path(the_file)
#
# extract the folder from the zipfile
#

with ZipFile(the_file, 'r') as myzip:
    myzip.extractall()
soundings_dir = Path('std_soundings')
metafile= soundings_dir / Path('metadata.json')
with open(metafile,'r') as meta_in:
    attribute_dict = json.load(meta_in)

filepaths = soundings_dir.glob('**/*.csv')
sound_dict={}
for the_path in filepaths:
    key = the_path.stem
    df = pd.read_csv(str(the_path))
    sound_dict[key] = df


# ## Data layout
# 
# There are five soundings.
# The soundings have six columns and 33 rows (i.e. 33 height levels).  The variables are
#  z, press, temp, rmix, den, o3den -- where rmix is the mixing ratio of water vapor, den is the dry air density and o3den is the ozone density.  The units are 
#  m, pa, K, kg/kg, kg/m^3, kg/m^3
#  
#  I will read the 6 column soundings into a [pandas (panel data) DataFrame](http://pandas.pydata.org/pandas-docs/stable/dsintro.html), which is like a matrix except the columns can be accessed by column name in addition to column number.  The main advantage for us is that it's easier to keep track of which variables we're plotting
#  

# In[3]:


list(sound_dict.keys())


# We use these keys to get a dataframe with 6 columns, and 33 levels.  Here's an example for the midsummer sounding

# In[4]:


midsummer=sound_dict['midsummer']
list(midsummer.columns)
midsummer.head()


# ### Plot  temp and vapor mixing ratio rmix ($\rho_{H2O}/\rho_{air}$)

# In[5]:


get_ipython().run_line_magic('matplotlib', 'inline')
meters2km=1.e-3
plt.close('all')
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(11,8))
for a_name,df in sound_dict.items():
    ax1.plot(df['temp'],df['z']*meters2km,label=a_name)
    ax1.set(ylim=(0,40),title='Temp soundings',ylabel='Height (km)',
       xlabel='Temperature (K)')

    ax2.plot(df['rmix']*1.e3,df['z']*meters2km,label=a_name)
    ax2.set(ylim=(0,8),title='Vapor soundings',ylabel='Height (km)',
       xlabel='vapor mixing ratio (g/kg)')
ax1.legend()
ax2.legend();


# ## Calculating scale heights for temperature and air density

# Here is equation 5 of the [hydrostatic balance notes](https://clouds.eos.ubc.ca/~phil/courses/atsc500/docs/hydro.pdf)

# $$\frac{ 1}{\overline{H_p}} =  \overline{ \left ( \frac{1 }{H} \right )} = \frac{\int_{0 }^{z}\!\frac{1}{H} dz^\prime  }{z-0} $$
# 
# where
# 
# $$H=R_d T/g$$

# and here is the Python code to do that integral:

# In[6]:


g=9.8  #don't worry about g(z) for this exercise
Rd=287.  #kg/m^3
def calcScaleHeight(T,p,z):
    """
    Calculate the pressure scale height H_p
    
    Parameters
    ----------
    
    T: vector (float)
      temperature (K)
      
    p: vector (float) of len(T)
      pressure (pa)
      
    z: vector (float) of len(T
      height (m)
      
    Returns
    -------
    
    Hbar: vector (float) of len(T)
      pressure scale height (m)
    
    """
    dz=np.diff(z)
    TLayer=(T[1:] + T[0:-1])/2.
    oneOverH=g/(Rd*TLayer)
    Zthick=z[-1] - z[0]
    oneOverHbar=np.sum(oneOverH*dz)/Zthick
    Hbar = 1/oneOverHbar
    return Hbar


# Similarly, equation (10) of the [hydrostatic balance notes](https://clouds.eos.ubc.ca/~phil/courses/atsc500/docs/hydro.pdf)
# is:
# 
# $$\frac{d\rho }{\rho} = - \left ( \frac{1 }{H} + 
#                    \frac{1 }{T} \frac{dT }{dz} \right ) dz \equiv - \frac{dz }{H_\rho} $$
#                    
# Which leads to 
# 
# $$\frac{ 1}{\overline{H_\rho}}  = \frac{\int_{0 }^{z}\!\left [ \frac{1}{H} + \frac{1 }{T} \frac{dT }{dz} \right ] dz^\prime  }{z-0} $$
# 
# and the following python function:
# 

# In[7]:


def calcDensHeight(T,p,z):
    """
    Calculate the density scale height H_rho
    
    Parameters
    ----------
    
    T: vector (float)
      temperature (K)
      
    p: vector (float) of len(T)
      pressure (pa)
      
    z: vector (float) of len(T
      height (m)
      
    Returns
    -------
    
    Hbar: vector (float) of len(T)
      density scale height (m)
    """
    dz=np.diff(z)
    TLayer=(T[1:] + T[0:-1])/2.
    dTdz=np.diff(T)/np.diff(z)
    oneOverH=g/(Rd*TLayer) + (1/TLayer*dTdz)
    Zthick=z[-1] - z[0]
    oneOverHbar=np.sum(oneOverH*dz)/Zthick
    Hbar = 1/oneOverHbar
    return Hbar


# ### How do $\overline{H_p}$ and $\overline{H_\rho}$ compare for the tropical sounding?

# In[8]:


sounding='tropics'
#
# grab the dataframe and get the sounding columns
#
df=sound_dict[sounding]
z=df['z'].values
Temp=df['temp'].values
press=df['press'].values


# In[9]:


#
# limit calculation to bottom 10 km
#
hit=z<10000.
zL,pressL,TempL=(z[hit],press[hit],Temp[hit])
rhoL=pressL/(Rd*TempL)
Hbar= calcScaleHeight(TempL,pressL,zL)
Hrho= calcDensHeight(TempL,pressL,zL)
print(f"pressure scale height for the {sounding} sounding is {Hbar*1.e-3:5.2f} km")
print(f"density scale height for the {sounding} is {Hrho*1.e-3:5.2f} km")


# ### How well do these average values represent the pressure and density profiles?

# In[10]:


theFig,theAx=plt.subplots(1,1)
theAx.semilogy(Temp,press/100.)
#
# need to flip the y axis since pressure decreases with height
#
theAx.invert_yaxis()
tickvals=[1000,800, 600, 400, 200, 100, 50,1]
theAx.set_yticks(tickvals)
majorFormatter = ticks.FormatStrFormatter('%d')
theAx.yaxis.set_major_formatter(majorFormatter)
theAx.set_yticklabels(tickvals)
theAx.set_ylim([1000.,50.])
theAx.set_title('{} temperature profile'.format(sounding))
theAx.set_xlabel('Temperature (K)')
theAx.set_ylabel('pressure (hPa)');


# Now check the hydrostatic approximation by plotting the pressure column against
# 
# $$p(z) = p_0 \exp \left (-z/\overline{H_p} \right )$$
# 
# vs. the actual sounding p(T):

# In[11]:


fig,theAx=plt.subplots(1,1)
hydroPress=pressL[0]*np.exp(-zL/Hbar)
theAx.plot(pressL/100.,zL/1000.,label='sounding')
theAx.plot(hydroPress/100.,zL/1000.,label='hydrostat approx')
theAx.set_title('height vs. pressure for tropics')
theAx.set_xlabel('pressure (hPa)')
theAx.set_ylabel('height (km)')
theAx.set_xlim([500,1000])
theAx.set_ylim([0,5])
tickVals=[500, 600, 700, 800, 900, 1000]
theAx.set_xticks(tickVals)
theAx.set_xticklabels(tickVals)
theAx.legend(loc='best');


# Again plot the hydrostatic approximation
# 
# $$\rho(z) = \rho_0 \exp \left (-z/\overline{H_\rho} \right )$$
# 
# vs. the actual sounding $\rho(z)$:

# In[12]:


fig,theAx=plt.subplots(1,1)
hydroDens=rhoL[0]*np.exp(-zL/Hrho)
theAx.plot(rhoL,zL/1000.,label='sounding')
theAx.plot(hydroDens,zL/1000.,label='hydrostat approx')
theAx.set_title('height vs. density for the tropics')
theAx.set_xlabel('density ($kg\,m^{-3}$)')
theAx.set_ylabel('height (km)')
theAx.set_ylim([0,5])
theAx.legend(loc='best');


# <a name="oct7assign"></a>
# 
# ### Assignment for Friday
# 
# Add cells to this notebook to:
# 
# 1\.  Print out the density and pressure scale heights for each of the five soundings
# 
# 2\.  Define a function that takes a sounding dataframe and returns the "total precipitable water", which is defined as:
# 
# $$W = \int_0^{z_{top}} \rho_v dz $$
# 
# Do a change of units to convert $kg\,m^{-2}$ to $cm\,m^{-2}$ using the density of liquid water (1000 $kg\,m^{-3}$) -- that is, turn the kg of water in the 1 square meter column into cubic meters and turn that into $cm/m^{-2}$
# 
# 3\.  Use your function to print out W for all five soundings
