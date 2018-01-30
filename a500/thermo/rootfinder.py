"""
   Two convenience functions for rootfinding
"""
import numpy as np
from scipy import optimize
from a405.utils.helper_funs import make_tuple

class BracketError(ValueError):
    """
    subclass ValueError so that I can add extra debugging
    information as the new member variable extra_info
    """
    def __init__(self,*args,**kwargs):
        super()  #call the base class constructor
        if 'extra_info' in kwargs:
            self.extra_info = kwargs['extra_info']

def find_interval(the_func, x, *args):
    """
    starting from a 2% difference, move out from a 
    point until the_func changes sign

    Parameters
    ----------

    the_func : function
               function that returns zero when on root
    
    x : float
        argument to the_func

    *args : tuple
            additional arguments for the_func

    Returns
    -------

    brackets : [left,right]
               left,right  brackets for root 
    """
    if x == 0.:
        dx = 1. / 50.
    else:
        dx = x / 50.

    maxiter = 40
    twosqrt = np.sqrt(2)

    failed = True
    for i in range(maxiter):
        dx = dx * twosqrt
        a = x - dx
        fa = the_func(a, *args)
        b = x + dx
        fb = the_func(b, *args)
        if (fa * fb < 0.):
            failed = False
            break
    if failed:
        #
        # load the debugging information into the BracketError exception as a
        # namedtuple
        #
        extra_info = make_tuple(dict(a=a,b=b,fa=fa,fb=fb,x=x,dx=dx,args=args))
        raise BracketError("Couldn't find a suitable range. Providing extra_info",extra_info=extra_info)
    return (a, b)


def fzero(the_func, root_bracket, *args, **parms):
    """
    simple wrapper for optimize.zeros.brenth

    Parameters
    ----------

    the_func : function
               function that returns zero when on root

    root_bracket : [left, right]
               left and right x values that bracket a sign change

    *args : tuple
            additional arguments for the_func

    **params: dict
              additional parameters for optimize.zeros.brenth

    Returns
    -------

    x value that produces the_func=0
    """
    answer = optimize.zeros.brenth(the_func,
                                   root_bracket[0],
                                   root_bracket[1],
                                   args=args,
                                   **parms)
    return answer


def test_rootfinder():
    """
    run unit tests for rootfinder
    """
    the_zero = fzero(np.sin, [12, 13]) * 180. / np.pi  #expecting 720 degrees
    np.testing.assert_almost_equal(the_zero, 720.)
    the_zero = fzero(np.sin, [18, 20], xtol=1.e-300, maxiter=80) * 180. / np.pi
    np.testing.assert_almost_equal(the_zero, 1080.)
    brackets = find_interval(np.sin, 25)
    the_zero = fzero(np.sin, brackets, xtol=1.e-300, maxiter=80) * 180. / np.pi
    np.testing.assert_almost_equal(the_zero, 1440.)


if __name__ == "__main__":
    test_rootfinder()
