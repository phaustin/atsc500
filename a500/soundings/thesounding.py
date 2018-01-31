import site
site.addsitedir('/Users/phil/repos/A405')
from importlib import reload
import pandas as pd

import readsoundings
reload(readsoundings)
sounding_dict=readsoundings.readsound(inFileName='soundings_july.txt')

with pd.HDFStore('out.h5','w') as store:
    for key,value in sounding_dict.items():
        thetime=key.strftime("Y%Y_%b_%d_%HZ")
        store.put(thetime,value,format='table')
        
with pd.HDFStore('out.h5','r') as store:
    print(store.keys())


with pd.HDFStore('out.h5','r') as store:
    sounding=store['Y2006_Jul_12_00Z']
    print(sounding)
