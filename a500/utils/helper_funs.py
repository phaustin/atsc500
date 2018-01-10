import numpy as np
from collections import namedtuple

def test_scalar(*args):
    """
    return true if every argument is a scalar

    parameters
    ==========

    args : list of arguments
           any number of comma-separated vector/scalar arguments

    Returns
    =======

    isscalar: bool
           true if every argument is type np.scalar
    """
    isscalar=True
    for item in args:
        isscalar = isscalar & np.isscalar(item)
    return isscalar

    
def make_tuple(in_dict,tupname='values'):
    """
    make a named tuple from a dictionary

    Parameters
    ==========

    in_dict: dictionary
         Any python object with key/value pairs

    tupname: string
         optional name for the new namedtuple type

    Returns
    =======

    the_tup: namedtuple
          named tuple with keys as attributes
    """
    the_tup = namedtuple(tupname, in_dict.keys())
    the_tup = the_tup(**in_dict)
    return the_tup

def find_centers(x):
    """
    return a vector of bin centers given the bin edges

    Parameters
    ==========

    x: numpy 1-d vector
       vector of edges of bins

    Returns
    =======
    center: numpy 1-d vector 
       vector of centers of bins
    
    """
    center = (x[1:] + x[:-1])/2.
    return center

