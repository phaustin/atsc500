import numpy as np
from add_path import add_path
add_path()


from a500.thermo import thermfuncs as tf
import pandas as pd

Temp=280. #K
press=90. #kPa
#
# find the saturation mixing ratio and potential temp from the temperature and pressure
#
print(tf.qs_tp(Temp,press))
print(tf.find_theta(Temp,press))
#
# find the dew point temperature and the lifting condensation level
#
psurf=100.  #kPa
Tsurf=290.
qvap=7.e-3  #kg/kg
Tdew = tf.tmr(qvap,psurf)
print('the surface dew point temp is: ',Tdew)
LCL = tf.LCL(Tdew,Tsurf,psurf)
print('the LCL is {:5.2f} kPa'.format(LCL))

#
# find thetal 
#
thetal = tf.alt_thetal(psurf,Tsurf,qvap)
#
# invert thetal for temperature, vapor, liquid
#
print(tf.t_uos_thetal(thetal,qvap,80.))
#
# make a sounding
#
press_levs=np.linspace(80,100.,20.)
press_levs=press_levs[::-1]
sounding=[]
for press in press_levs:
    sounding.append(tf.t_uos_thetal(thetal,qvap,press))
    
df_sounding=pd.DataFrame.from_records(sounding)

