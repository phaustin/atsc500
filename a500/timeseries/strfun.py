from a500.utils.data_read import download
import numpy as np
from matplotlib import pyplot as plt

import os
from numba import jit
import numba

@jit(nopython=True)
def calc_sf2(vector,num_seps):
    totvals=len(vector)
    accum=totvals
    halfvals = np.int32(totvals/2.)
    spacing = np.arange(1,halfvals,1,np.int32)
    accum= np.zeros(spacing.shape,dtype=np.float32)
    count = np.zeros(spacing.shape,np.int32)
    for accum_index,the_sep in enumerate(spacing[:num_seps]):
        vals=np.arange(the_sep,halfvals,1,np.int32)
        for vec_index in vals:
             accum[accum_index] = accum[accum_index] + (vector[vec_index] - vector[vec_index - the_sep])**2.
             count[accum_index] = count[accum_index] + 1
        accum[accum_index] = accum[accum_index]/count[accum_index]
    return accum,count


download('aircraft.npz',root='http://clouds.eos.ubc.ca/~phil/docs/atsc500/data')
data = np.load('aircraft.npz')
maxlength=int(1.e5)
vector=data['wvel'][:maxlength]
out=calc_sf2(vector,1000)    
        
        
