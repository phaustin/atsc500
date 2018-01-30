"""
  Module for functions to calculate vapor pressure,
  thetae, etc.
"""

from collections import Iterable
import numpy as np
import doctest

from a405.thermo.constants import constants as c
from a405.thermo import rootfinder as rf
import numpy.testing as ntest
from a405.thermo.rootfinder import BracketError
from .constants import constants as c

def find_lv(temp):
    """
    Calculates the temperature dependent
    enthalpy of evaporation

    Parameters
    ----------

    temp : float or array_like
           Temperature of parcel (K).

    Returns
    -------

    lv : float or list
         entalpy of evaporation


    Examples
    --------

    >>> find_lv(293.)
    2547052.0

    References
    ----------

    Day 13 moist static energy notes
    """
    lv = c.lv0 - (c.cpv - c.cl) * (temp - c.Tc)
    return lv


def find_esat(temp):
    """
    Calculates the saturation water vapor pressure over a flat
    surface of water at temperature 'temp'.

    Parameters
    ----------

    temp : float or array_like
           Temperature of parcel (K).

    Returns
    -------

    esatOut : float or list
        Saturation water vapour pressure (Pa).

    Examples
    --------

    >>> find_esat(300.)
    3534.5196668891358
    >>> find_esat([300., 310.])
    array([ 3534.5197,  6235.5322])

    References
    ----------
    Emanuel 4.4.14 p. 117
      
    """
    # determine if temp has been input as a vector
    is_scalar = True
    if isinstance(temp, Iterable):
        is_scalar = False
    temp = np.atleast_1d(temp)
    Tc = temp - 273.15
    esatOut = 611.2 * np.exp(17.67 * Tc / (Tc + 243.5))
    # if temp is a vector
    if is_scalar:
        esatOut = esatOut[0]
    return esatOut


def find_resid_rsat(Tguess, rsat, press):
    """

   Calculate residual between target rsat and guess to
   rootfind saturation temperature for constant rsat   
   
   Parameters
   ----------
      
   Tguess: float
           (K) guess temperature from rootfinder 

   rsat:   float
           (kg/kg) target saturation mixing ratio
            
   press:  float
           (hPa) pressure level for rootfind

   Returns
   ------- 
   
   residual: float
             (kg/kg) differnce between target and guess rsat

   Reference
   ---------
      
   see thompkins 2.20
     
    """
    esat = find_esat(Tguess) * 0.01  #convert to hPa
    residual = rsat - c.eps * esat / (press - esat)
    return residual


def find_rsat(temp, press):
    """
  
   calculate the saturation mixing ratio (kg/kg) at (temp,press)

   Parameters
   --------- 
       
   temp: float
         temperature (K) 

   press: float
         pressure (Pa)

   Returns
   -------
       
    rsat: float
          satruation mixing ratio  (kg/kg)
    """
    esat = find_esat(temp)
    rsat = c.eps * esat / (press - esat)
    return rsat


def tinvert_rsat(Tstart, rsat, press):
    """
    rootfind the temp that produces rsat at press.

    Parameters
    ----------

    temp : float
           temperature (K)

    rsat : float
           saturation mixing ratio (kg/kg)

    press : float 
            pressure (hPa)

    Returns
    -------

    Tdew : temperature (K) at with air is saaturated
    """
    brackets = rf.find_interval(find_resid_rsat, Tstart, rsat, press)
    temp = rf.fzero(find_resid_rsat, brackets, rsat, press)
    return temp


def find_theta(temp, press, rv=0):
    """
    Computes potential temperature.
    Allows for either temp,p or T,p,rv as inputs.

    Parameters
    ----------

    temp : float
        Temperature (K).

    press : float
        Pressure (Pa).

    rv : float, optional
        Vapour mixing ratio (kg,kg). Can be appended as an argument
        in order to increase precision of returned 'theta' value.


    Returns
    -------

    thetaOut : float
        Potential temperature (K).

    
    Raises
    ------

    NameError
        If an incorrect number of arguments is provided.
    
    
    References
    ----------
    Emanuel p. 111 4.2.11


    Examples
    --------
    >>> find_theta(300., 8.e4) # Only 'temp' and 'press' are input.
    319.72798180767984
    >>> find_theta(300., 8.e4, rv=0.001) # 'temp', 'press', and 'rv' all input.
    319.72309475657323
    
    """

    power = c.Rd / c.cpd * (1. - 0.24 * rv)
    thetaOut = temp * (c.p0 / press)**power
    return thetaOut


def convertTempToSkew(Temp, press, skew):
    """
    convertTempToSkew(Temp, press, skew)

    Determines the transformed temperature in plotting coordinates.
    
    Parameters
    ----------
    Temp : float
        Temperature (degC)
    press : float
        Pressure (hPa).
    skew : int
        Designated skew factor of temperature.

    Returns
    ----
    tempOut : float
        Converted temperature (degC).

    Examples
    --------
    >>> convertTempToSkew(30., 800. , 30)
    -170.53835183003781
    
    """

    tempOut = Temp - skew * np.log(press)
    return tempOut


def convertSkewToTemp(xcoord, press, skew):
    """
    convertSkewToTemp(xcoord, press, skew)

    Determines temperature from knowledge of a plotting coordinate
    system and corresponding plot skew.
    
    Parameters
    ----------
    xcoord : int
        X coordinate in temperature plotting coordinates.
    press : float
        Pressure (hPa).
    skew : int
        Skew of a given coordinate system.

    Returns
    ----
    Temp : float
        Converted temperature in degC.

    Examples
    --------
    >>> convertSkewToTemp(-170.53835183003781, 800., 30)
    30.0
    
    """
    Temp = xcoord + skew * np.log(press)
    return Temp


def find_thetaes(Temp, press):
    """
    Calculates the true equivalent potential temperature of an air
    parcel assuming saturation

    Parameters
    ----------
    Temp : float
        Temperature (K).
    press : float
        Pressure (Pa).

    Returns
    ----
    thetaep : float
        Pseudo equivalent potential temperature (K).


    Notes
    -----
    Empirical fit to Emanuel 4.7.8 -- neglects heat
    capacity of liquid water and ice


    References
    ----------
    Emanuel 4.7.9 p. 132


    Examples
    --------

    >>> find_thetaes(300., 8.e4)
    399.53931578042267
    """
    # The parcel is saturated - prohibit supersaturation with Td > T.
    Td = Temp
    rv = find_rsat(Td, press)
    thetaep = find_thetaet(Td, rv, Temp, press)
    #
    # peg this at 450 so rootfinder won't blow up
    #
    if thetaep > 450.:
        thetaep = 450
    return thetaep


def find_Tv(temp, rvap, rl=0.):
    """
    Calculate the density (virtual) temperature


    Parameters
    ----------

    temp : float or np.array
        temperatue (K).

    rvap: float or np.array
        vapor  mixing ratio (kg/kg)

    rl : float or np.array
        liquid mixing ratio (kg/kg)

    Returns
    -------

    Tv : float or np.array
        density temperature (K).

    References
    ----------
    Thompkins eq. 2.38


    Examples
    --------
    >>> find_Tv(300.,1.e-2)# Parcel is unsaturated.
    301.866
    >>> find_Tv(280.,1.e-2, 1.e-3)  #Parcel is saturated
    281.4616
    """
    return temp * (1 + c.eps * rvap - rl)


def find_thetaet(Td, rt, T, p):
    """
    Calculates the true equivalent potential temperature of a
    parcel from its entropy


    Parameters
    ----------

    Td : float
        Dewpoint temperature (K).

    rt: float
        total water mixing ratio (kg/kg)

    T : float
        Temperature (K).

    p : float
        Pressure (Pa).


    Returns
    -------

    thetaetOut : float
        true equivalent potential temperature (K).


    References
    ----------
    Emanuel 4.5.11 p. 120 or my Day 10 equivalent potential temperature notes


    Examples
    --------
    >>> find_thetaet(300.,1.e-2, 280., 8.e4) # Parcel is saturated.
    319.3543617606212
    >>> find_thetaet(280.,1.e-2, 300., 8.e4)  #Parcel is unsaturated
    342.52792353970784
    """
    if Td < T:
        #parcel is unsaturated
        [Tlcl, plcl] = find_lcl(Td, T, p)
        rv = find_rsat(Td, p)
    else:
        #parcel is saturated -- prohibit supersaturation with Td > T
        Td = T
        rv = find_rsat(T, p)
    e = find_esat(Td)
    esat = find_esat(T)
    vapor_term = rv * c.Rv * np.log(e / esat)
    #
    # turn off water vapor if not in liquid water saturation
    # domain
    #
    if np.isinf(vapor_term) or (e > 1.e5) or (esat > 1.e5) \
       or (e < 1) or (esat < 1):
        vapor_term = 0
        e = 0
    #print('thetaet: ',e,esat,T,Td,p,vapor_term)
    pd = p - e  # dry
    cp = c.cpd + rt * c.cl
    lv = find_lv(T)
    s = cp * np.log(T) - c.Rd * np.log(pd) + lv * rv / T - vapor_term
    logthetae = (s + c.Rd * np.log(c.p0)) / cp
    thetaet = np.exp(logthetae)
    #
    # peg this at 450 so rootfinder won't blow up
    #
    if (thetaet > 450.):
        thetaet = 450
    if (thetaet < 100.):
        thetaet = 100.
    return thetaet


def find_thetaep(Td, T, p):
    """
    Calculates the pseudo equivalent potential temperature of a
    parcel using Bolton's formula


    Parameters
    ----------

    Td : float
        Dewpoint temperature (K).

    T : float
        Temperature (K).

    p : float
        Pressure (Pa).


    Returns
    -------

    thetaepOut : float
        Pseudo equivalent potential temperature (K).


    Notes
    -----
    Note that the pseudo equivalent potential temperature of an air
    parcel is not a conserved variable.


    References
    ----------
    Emanuel 4.7.9 p. 132


    Examples
    --------
    # note difference between true thetae (find_thetaes and find_thetaet)
    # and pseudo-adiabat (find_thetaep)

    >>> find_thetaep(280., 300., 8.e4) # Parcel is unsaturated.
    344.99830738253371

    >>> rt = find_rsat(280, 8.e4)
    >>> find_thetaet(280., rt, 300., 8.e4) # Parcel is unsaturated.
    342.93068265625539


    >>> find_thetaep(280., 280., 8.e4) # Parcel is saturated.
    321.5302927767795


    >>> find_thetaes(280., 8.e4) # Parcel is saturated.
    319.72687107952828

    """
    if Td < T:
        #parcel is unsaturated
        [Tlcl, plcl] = find_lcl(Td, T, p)
        rv = find_rsat(Td, p)
    else:
        #parcel is saturated -- prohibit supersaturation with Td > T
        Tlcl = T
        rv = find_rsat(T, p)

    thetaval = find_theta(T, p, rv)
    thetaepOut = thetaval * np.exp(rv * (1 + 0.81 * rv) \
                                   * (3376. / Tlcl - 2.54))
    #
    # peg this at 450 so rootfinder won't blow up
    #

    if (thetaepOut > 450.):
        thetaepOut = 450

    return thetaepOut


def find_lcl(Td, T, p):
    """
    find_lcl(Td, T, p)

    Finds the temperature and pressure at the lifting condensation
    level (LCL) of an air parcel.

    Parameters
    ----------
    Td : float
        Dewpoint temperature (K).

    T : float
        Temperature (K).

    p : float
        Pressure (Pa)

    Returns
    -------

    Tlcl : float
        Temperature at the LCL (K).
    plcl : float
        Pressure at the LCL (Pa).

    Raises
    ------

    NameError
        If the air is saturated at a given Td and T (ie. Td >= T)
    
    Examples
    --------

    >>> [Tlcl, plcl] =  find_lcl(280., 300., 8.e4)
    >>> print(np.array([Tlcl, plcl]))
    [   275.7625  59518.9287]
    >>> find_lcl(300., 280., 8.e4)
    Traceback (most recent call last):
        ...
    NameError: parcel is saturated at this pressure

    References
    ----------
    Emanuel 4.6.24 p. 130 and 4.6.22 p. 129
    
    """
    hit = Td >= T
    if hit is True:
        raise NameError('parcel is saturated at this pressure')

    e = find_esat(Td)
    ehPa = e * 0.01
    #Bolton's formula requires hPa.
    # This is is an empircal fit from for LCL temp from Bolton, 1980 MWR.
    Tlcl = (2840. / (3.5 * np.log(T) - np.log(ehPa) - 4.805)) + 55.

    r = c.eps * e / (p - e)
    cp = c.cpd + r * c.cpv
    logplcl = np.log(p) + cp / (c.Rd * (1 + r / c.eps)) * \
              np.log(Tlcl / T)
    plcl = np.exp(logplcl)

    return Tlcl, plcl


def find_rvrl(Temp, rT, press):
    """
    Computes the vapour and liquid water mixing ratios.

    Parameters
    ----------
    Temp : float
        Temperature (K).
    rT : float
        Total water mixing ratio (kg/kg).
    press : float
        Pressure (Pa).


    Returns
    ----
    rv : float
        Water vapour mixing ratio (kg/kg).
    rl : float
        Liquid water mixing ratio (kg/kg).


    Raises
    ----
    AssertionError
        If any of the inputs are in vector form.

    Examples
    --------

    >>> print(np.array(find_rvrl(250., 0.01, 8.e4))*1.e3)
    [ 0.7433  9.2567]

    >>> print(np.array(find_rvrl(305., 0.01, 8.e4))*1.e3)
    [ 10.   0.]

    """
    rsVal = find_rsat(Temp, press)
    if rsVal > rT:  #unsaturated
        rv = rT
        rl = 0
    else:  #saturated
        rv = rsVal
        rl = rT - rv
    return rv, rl


def find_Tmoist(thetaE0, press):
    """
    Calculates the temperatures along a moist adiabat.
    
    Parameters
    ----------

    thetaE0 : float
        Initial equivalent potential temperature (K).
    press : float or array_like
        Pressure (Pa).

    Returns
    -------
    Temp : float or array_like
        Temperature (K) of thetaE0 adiabat at 'press'.

    Examples
    --------
    >>> np.array([find_Tmoist(300., 8.e4)])
    array([ 271.0638])
    >>> find_Tmoist(330., 8.e4)
    283.7226584032411
    """
    Tstart = c.Tc
    try:
        brackets = rf.find_interval(thetaes_diff, Tstart, thetaE0, press)
        Temp = rf.fzero(thetaes_diff, brackets, thetaE0, press)
    except BracketError as e:
        print("couldn't find bracket: debug info: ", e.extra_info)
        Temp = np.nan
    return Temp


def thetaes_diff(Tguess, thetaE0, press):
    """
    use true thetae (thetaes) for rootfinder

    Parameters
    ----------
    Tguess : float
        Trial temperature value (K).
    ws0 : float
        Initial saturated mixing ratio (kg/kg).
    press : float
        Pressure (Pa).

    Returns
    -------
    theDiff : float
        The difference between the values of 'thetaEguess' and
        'thetaE0'. This difference is then compared to the tolerance
        allowed by brenth.
        
    """
    thetaes_guess = find_thetaes(Tguess, press)

    #when this result is small enough we're done
    the_diff = thetaes_guess - thetaE0
    return the_diff


def thetaep_diff(Tguess, thetaE0, press):
    """
    use pseudo thetae (thetaep) for rootfinder

    Parameters
    ----------
    Tguess : float
        Trial temperature value (K).
    ws0 : float
        Initial saturated mixing ratio (kg/kg).
    press : float
        Pressure (Pa).

    Returns
    -------
    theDiff : float
        The difference between the values of 'thetaEguess' and
        'thetaE0'. This difference is then compared to the tolerance
        allowed by brenth.
        
    """
    thetaes_guess = find_thetaep(Tguess, Tguess, press)

    #when this result is small enough we're done
    the_diff = thetaes_guess - thetaE0
    return the_diff


def tinvert_thetae(thetaeVal, rT, press):
    """
    temp,rv,rl=tinvert_thetae(thetaeVal, rT, press)

    Uses a rootfinder to determine the temperature for which the
    pseudo equivilant potential temperature (thetaepress) is equal to the
    equivilant potential temperature (thetae) of the parcel.

    Parameters
    ----------
    thetaeVal : float
        Thetae of parcel (K).
    rT : float
        Total water mixing ratio (kg/kg).
    press : float
        Pressure of parcel in (Pa).

    Returns
    -------

    theTemp : float
        Temperature for which thetaep equals the parcel thetae (K).
    rv : float
        Vapor mixing ratio of the parcel (kg/kg).
    rl : float
        liquid water mixing ratio of the parcel (kg/kg) at 'press'.

    Raises
    ------
    IOError
        If 'press' is larger than 100000 Pa.

    Examples
    --------

    >>> tinvert_thetae(300., 0.001, 8.e4)
    (278.683729619619, 0.001, 0)
    """
    if press > 1.e5:
        raise IOError('expecting pressure level less than 100000 Pa')
    # The temperature has to be somewhere between thetae
    # (T at surface) and -40 deg. C (no ice).
    Tstart = c.Tc
    brackets = rf.find_interval(find_resid_thetae, Tstart, thetaeVal, rT,
                                press)
    theTemp = rf.fzero(find_resid_thetae, brackets, thetaeVal, rT, press)
    rv, rl = find_rvrl(theTemp, rT, press)
    return theTemp, rv, rl


def find_resid_thetae(Tguess, thetaeVal, rT, press):
    """
   Calculate residual between target thetaeVal and guess to
   rootfind temperature for constant thetae
   
   Parameters
   ----------
      
   Tguess: float
           (K) guess temperature from rootfinder 

   rsat:   float
           (kg/kg) target saturation mixing ratio
            
   press:  float
           (hPa) pressure level for rootfind

   Returns
   ------- 
   
   residual: float
             (kg/kg) differnce between target and guess rsat

   References
   ----------
      
   see thompkins 2.20
   """

    rv, rl = find_rvrl(Tguess, rT, press)
    tdGuess = find_Td(rv, press)
    # Iterate on Tguess until this function is
    # zero to within tolerance.
    return thetaeVal - find_thetaet(tdGuess, rT, Tguess, press)


def find_buoy(adia_Tv, env_Tv):
    """
    Calculates the buoyancy given an parcel and environment virtural temp

    Parameters
    ----------

    adia_Tv : float
        density temperature of the parcel (K)
    env_Tv : float
        density temperature of environment (K)

    Returns
    -------

    buoy : float
        buoyancy (m/s/s)

    Examples
    --------

    >>> find_buoy(286.,285.)
    0.03438596491228071

    References
    ----------

    Thompkins equation 3.3
    
    """
    buoy = c.g0 * (adia_Tv - env_Tv) / env_Tv
    return buoy


def find_thetal(press, temp, rt):
    """
    Calculates the true equivalent potential temperature of an air
    parcel assuming saturation

    Parameters
    ----------
    temp : float
        Temperature (K).
    press : float
        Pressure (Pa).

    rt : float
       mixing ratio  (kg/kg

    Returns
    ----
    thetal : float
        liquid water potential temperature (K).


    Notes
    -----
    Empirical fit 


    References
    ----------
    Emanuel 4.5.15 p. 121

    """
    rsat = find_rsat(temp, press)
    saturated = rt > rsat
    if saturated:
        rl = rt - rsat
    else:
        rl = 0.
    lv = find_lv(temp)
    theta = find_theta(temp, press)
    thetal = theta * np.exp(-lv * rl / (c.cpd * temp))
    return thetal


def find_Td(rv, press):
    """
    Calculates the due point temperature of an air parcel.

    Parameters
    ----------

    rv : float
        Mixing ratio (kg/kg).
    press : float
        Pressure (Pa).

    Returns
    -------

    Td : float
        Dew point temperature (K).

    Examples
    --------

    >>> find_Td(0.001, 8.e4)
    253.39429263963504

    References
    ----------

    Emanuel 4.4.14 p. 117
    
    """
    e = rv * press / (c.eps + rv)
    denom = (17.67 / np.log(e / 611.2)) - 1.
    Td = 243.5 / denom
    Td = Td + 273.15
    return Td


def test_therm():
    """
   execute unit tests for thermlib
   """
    Tlcl, plcl = find_lcl(280., 300., 8.e4)
    ntest.assert_almost_equal(Tlcl, 275.7625, decimal=3)
    ntest.assert_almost_equal(plcl, 59518.928, decimal=2)
    ntest.assert_almost_equal(
        find_thetaep(280., 300., 8.e4),
        344.998307,
        decimal=5)  # Parcel is unsaturated.
    ntest.assert_almost_equal(
        find_thetaep(300., 280., 8.e4),
        321.53029,
        decimal=5)  # Parcel is saturated.
    ntest.assert_almost_equal(find_esat(300.), 3534.51966, decimal=2)
    ntest.assert_almost_equal(find_thetaes(300., 8.e4), 399.53931, decimal=4)
    ntest.assert_allclose(find_esat([300., 310.]), [3534.51966, 6235.53218])
    ntest.assert_almost_equal(find_Tmoist(300., 8.e4), 271.063785, decimal=4)
    ntest.assert_almost_equal(find_Tmoist(330., 8.e4), 283.722658, decimal=4)
    ntest.assert_almost_equal(find_Tv(300., 1.e-2), 301.866, decimal=3)
    ntest.assert_almost_equal(find_Tv(280., 1.e-2, 1.e-3), 281.4616, decimal=3)


if __name__ == "__main__":
    np.set_printoptions(precision=4)
    test_therm()
    doctest.testmod()
