"""
In this library, I try to use the least approximation in the
thermodynamic calculation. Some necessary approximation is put at the
beginning. For functions with input parameters less than 3, the
function name is self-explaining; the string before the underline is
the output variable while the string after the underline indicates the
input variables. The convention is t: temperature, tv: virtual
temperature, tro: density temperature, th: potential temperature, thv:
virtual potential temperature, thro:virtual density temperature, p:
pressure, r: water vapor mixing ratio, rt: total water mixing ratio,
q: specific humidity, qt: specific total water content, ql: specific
liquid water content, e: water vapor pressure, es: saturated water
vapor pressure. For functions with input parameters more than 3
parameters, please look at the detail decription. For functions doing
similar calculation but with different approximations the second
underline is given and the string after the second underline indicates
the approximation name.

"""
import numpy as np
from . import thermconst as tc
from .thermconst import LV0,CPD
import textwrap
import sys
from ..utils.helper_funs import test_scalar

def alt_thetal(press,Temp,qt):
    press=np.atleast_1d(press)
    Temp=np.atleast_1d(Temp)
    qt=np.atleast_1d(qt)
    qsat=qs_tp(Temp,press)
    ql=qt - qsat
    ql[ql < 0]=0.
    theta = Temp*(100./press)**0.288
    L=2.5e6
    cpd=1004.
    thetal=theta - (theta/Temp)*(L/cpd)*ql
    return thetal

def thetal(press,Temp,rt):
    """
    input: press in kPa  (1d vector or scalar)
           temp in K  (1d vector or scalar
           mixing ratio in kg/kg  (1d vector or scalar)
    # emanuel p. 121  4.5.15
    """
    isscalar = test_scalar(press,Temp,rt)
    press=np.atleast_1d(press)
    Temp=np.atleast_1d(Temp)
    rt = np.atleast_1d(rt)
    rl=np.empty_like(rt)
    rsat=rs_tp(Temp,press)
    CPn=tc.RD + rt*tc.RV
    chi=CPn/(tc.CPD + rt*tc.CPV)
    gamma=(rt*tc.RV)/CPn
    saturated=rt > rsat
    rl[saturated]=rt[saturated] - rsat[saturated]
    unsaturated=np.logical_not(saturated)
    rl[unsaturated]=0
    LV=L_t(Temp)
    theta=Temp*(tc.P0/press)**chi
    term1=(1. - rl/(tc.EPS + rt))**chi
    term2=(1 - rl/rt)**(-gamma)
    term3= -LV*rl/(CPn*Temp)
    term3=np.exp(term3)
    theThetal=theta*term1*term2*term3
    if isscalar:
        theThetal = theThetal[0]
    return theThetal

def find_theta(Temp, press):
    """
    Parameters
    ----------
    
    Temp: float
      temperature (K)

    press: float
      pressure (kPa)

    Returns
    -------

    theta: float
       potential temperature (K)
    """
    theta = Temp*(100./press)**0.288
    return theta


def es_t_boton(Temp):
    """in: r=temperature in Kelvin
       out:es=saturation vapour pressure over water (kPa)
    """
    ttmp=Temp-273.15
    result = 0.1*6.112*np.exp(17.67*ttmp/(ttmp+243.5))
    return result

def L_t(T):
    """in: T=temperature in Kelvin
       out:L=temperature dependent latent heat of vaporization
    """
    L=tc.LV0 - tc.CPVMCL*(T - tc.T0)
    return L

def q_r(r):
    """in: r=mixing ratio r (kg/kg)
       out:q=specific humidity (kg/kg)
    """
    result = r/(1+r)
    return result

def r_q(q):
    """in: q=specific humidity (kg/kg)
       out:r=mixing ratio (kg/kg)
    """
    result = q/(1-q)
    return result

def rs_tp(Temp,press):
    """in: Temp=temperature (K),press=pressure (kPa)
       out:rs=saturation mixing ratio (kg/kg)
    """
    theEs=es_t_boton(Temp)
    result = tc.EPS*theEs/(press-theEs);
    return result

def dqs_dt(temp,p):
   """in: temp=temperature (K),p= pressure (kPa)
      out: derivive of the vapor pressure wrt temperature
           (kPa/K)
   """
   lv=L_t(temp)
   qs=qs_tp(temp,p)
   dqsdt = qs*lv/(tc.Rv*temp**2)
   return dqsdt
   
    
      

def qs_tp(t,p):
    """in: t=temperature (K),p= pressure (kPa)
       out:rs=saturation specific humidity (kg/kg)
    """
    rs=rs_tp(t,p)
    result =q_r(rs)
    return result

def t_thpr(th,p,r):
    """in: th=potential temperature (K),p= pressure (kPa)
           r=water vapor mixing ratio (kg/kg)
       out:t=temperature (K)
    """
    th=np.asarray(th)
    p=np.asarray(p)
    r=np.asarray(r)
    Rp=tc.RD*(1.+r/tc.EPS)/(1.+r)
    CPp=tc.CPD*(1.+r*(tc.CPV/tc.CPD))/(1.+r)
    result=th * (p/tc.P0)**(Rp/CPp)
    return result

def tv_tr(t,r):
    """in: t=temperature (K),r= water vapor mixing ration (kg/kg)
       out:rv=virtual potential temperature (K)
    """
    result  = t * ( 1.+ r / tc.EPS ) / ( 1.+ r )
    return result

def thv_tvp(tv,p):
    """in: tv=virtual temperature (K),p= pressure (kPa)
       out:rs=virtual potential temperature (kg/kg)
    """
    result  = tv*(tc.P0/p)**tc.RDOCP
    return result

def tro_trtp(t,rt,pressKpa):
    """in: t=temperature (K),pressKpa= pressure (kPa)
           rt=mixing ratio (kg/kg)
       out:tro=density temperature (K)  (vector)
    """
    isscalar = test_scalar(t,rt,pressKpa)
    if not isscalar:
        t=np.atleast_1d(t)
        rt=np.atleast_1d(rt)
        rv=np.empty_like(rt)
        pressKpa=np.atleast_1d(pressKpa)
    rs=rs_tp(t,pressKpa)
    out=np.empty(rs.shape,dtype=np.float)
    saturated = (rt>=rs)
    unsaturated=np.logical_not(saturated)
    rl=np.zeros(out.shape,dtype=np.float)
    rl[saturated]=rt[saturated] - rs[saturated]
    rl[unsaturated]=0
    rv[saturated]=rs[saturated]
    rv[unsaturated]=rt[unsaturated]
    out[saturated]=t[saturated]*(1.+0.611*rs[saturated]-rl[saturated])
    out[unsaturated]=t[unsaturated]*(1.+ 0.611*rt[unsaturated])
    out=t*(1 + rv/tc.EPS)/(1. + rt)
    if isscalar:
        out=out[0]
    return out

from scipy import optimize

def findDiffT(Tguess,htarget,rtTarget,press,gzval):
    rsat=rs_tp(Tguess,press)
    CPn=(tc.CPD + rtTarget*tc.CPV)
    gzvalp=(1+rtTarget)*gzval
    #CPn=tc.CPD
    rl=rtTarget - rsat
    if rl < 0:
        rl=0.
    if rtTarget > rsat:
        LV=L_t(Tguess)
        hguess=CPn*Tguess - LV*rl + gzvalp
        sat=True
    else:
        hguess=CPn*Tguess + gzvalp
        sat=False
    #print "in findDiffT: ",Tguess,htarget,rtTarget,p,gz,rs,hguess,sat
    theDiff=htarget - hguess
    return theDiff

def findDiffTthetal(Tguess,thetaltarget,rtTarget,press):
    thetalguess=thetal(press,Tguess,rtTarget)
    #print "in findDiffT: ",Tguess,htarget,rtTarget,p,gz,rs,hguess,sat
    theDiff=thetaltarget - thetalguess
    return theDiff



def t_uos(h,rt,p,gz):
    """in: h=liquid water static energy per unit of moist air h(J/kg)
           rt=total water mixing ratio (kg/kg)
           p=pressure (kPa),gz=geopotential height gz (m2/s2)
       out:dictionary with t, ql and x
       Emanuel p. 123, 4.5.25
    """
    ## def findDiffT(Tguess,htarget,qtTarget,p,gz):
    ##     rs=rs_tp(Tguess,p)
    ##     if qtTarget > rs:
    ##         CPn=(tc.CPD + qtTarget*tc.CL)
    ##         CPn=tc.CPD
    ##         LV=L_t(Tguess)
    ##         hguess=CPn*Tguess -LV*(qtTarget- rs) + gz
    ##         sat=True
    ##     else:
    ##         CPn=(tc.CPD + qtTarget*tc.CL)
    ##         CPn=tc.CPD
    ##         hguess=CPn*Tguess + gz
    ##         sat=False
    ##     #print "in findDiffT: ",Tguess,htarget,qtTarget,p,gz,rs,hguess,sat
    ##     theDiff=htarget - hguess
    ##     return theDiff



    Tlow=230.
    Thigh=340.
    t1=optimize.zeros.brenth(findDiffT, Tlow, Thigh, (h,rt, p,gz));
    result={}
    rs1=rs_tp(t1,p); 
    if   (rt>rs1):
        result["T"]=t1; result["RL"]=rt-rs1; result["X"]=1;result["RV"]=rs1
    else:
        result["T"]=t1; result["RL"]=0.0;    result["X"]=0;result["RV"]=rt
    return result


def t_uos_thetal(thetal,rt,p):
    """in: h=liquid water static energy per unit of moist air h(J/kg)
           rt=total water mixing ratio (kg/kg)
           p=pressure (kPa),gz=geopotential height gz (m2/s2)
       out:dictionary with t, ql and x
       Emanuel p. 123, 4.5.25
    """
    Tlow=230.
    Thigh=340.
    t1=optimize.zeros.brenth(findDiffTthetal, Tlow, Thigh, (thetal,rt, p));
    result={}
    rs1=rs_tp(t1,p)
    if   (rt>rs1):
        result["T"]=t1; result["RL"]=rt-rs1; result["X"]=1;result["RV"]=rs1
    else:
        result["T"]=t1; result["RL"]=0.0;    result["X"]=0;result["RV"]=rt
    return result

def sat_line(h,p,gz):
    """find rtList that is exactly saturated at static energy
       values given by hlist
    """

    def findDiffR(rtGuess,htarget,p,gz):
        #parcel is just saturated
        rv=rtGuess
        Tguess=tmr(rtGuess,p)
        CPn=tc.CPD + tc.CPV*rv
        gzp=gz*(1 + rv)
        hguess=CPn*Tguess + gzp
        theDiff=htarget - hguess
        return theDiff

    rlow=1.e-4
    rhigh=85.e-3
    rsat=optimize.zeros.brenth(findDiffR, rlow, rhigh, (h,p,gz),xtol=1.e-7);
    return rsat


def sat_line_thetal(thetalTarget,p):
    """find rtList that is exactly saturated at static energy
       values given by hlist
    """

    def findDiffR(rtGuess,thetalTarget,p):
        #parcel is just saturated
        rv=rtGuess
        Tguess=tmr(rtGuess,p)
        thetalGuess=thetal(p,Tguess,rtGuess)
        theDiff=thetalTarget - thetalGuess
        #print "debug: ",thetalTarget,rtGuess,theDiff
        return theDiff

    rlow=1.e-4
    rhigh=85.e-3
    rsat=optimize.zeros.brenth(findDiffR, rlow, rhigh, (thetalTarget,p),xtol=1.e-7);
    return rsat



def calc_dens(h,rt,p,gz):
    result=t_uos(h,rt,p,gz)
    T=result["T"]
    rl=result["RL"]
    rv=result["RV"]
    #print "calc_dens: here is h,rt,T,rl,qv: ",h,rt,T,rl,qv
    #4.3.5 p. 113
    Dguess=p*1.e3*(rl + rv)/(tc.RD*T*(1. + (rv/tc.EPS)))
    return Dguess

def calc_dens_thetal(thetal,rt,p):
    result=t_uos_thetal(thetal,rt,p)
    T=result["T"]
    rl=result["RL"]
    rv=result["RV"]
    #print "calc_dens: here is h,rt,T,rl,qv: ",h,rt,T,rl,qv
    #4.3.5 p. 113
    Dguess=p*1.e3*(rl + rv)/(tc.RD*T*(1. + (rv/tc.EPS)))
    return Dguess


output="""
    r %(r)f
    T %(T)f
    Dtarget %(Dtarget)f
    hmGuess %(hmGuess)f
    p %(p)f
    gz %(gz)f
    Dguess %(Dguess)f
    theDiff %(theDiff)f
    """

output=textwrap.dedent(output)

def findDiffD(hmGuess,Dtarget,r,p,gz):
    #parcel is just saturated
    #print "in diffD 1, ",hmGuess,Dtarget
    result=t_uos(hmGuess,r,p,gz)
    #print "in diffD 2, ",hmGuess,Dtarget
    T=result["T"]
    rl=result["RL"]
    rv=result["RV"]
    #print "calc_dens: here is h,rt,T,rl,qv: ",h,rt,T,rl,qv
    #4.3.5 p. 113
    Dguess=p*1.e3*(1. + rl + rv)/(tc.RD*T*(1. + (rv/tc.EPS)))
    #print "checkit"
    theDiff=Dtarget - Dguess
    #print "in diffD, ",hmGuess,result,theDiff
    ## if result['X']==0 and theDiff < 1.e-5:
    ##     print "in findDiff, unsat: "
    ##     print output % {'r':r,'T':T,'Dtarget':Dtarget,'p':p,'gz':gz,
    ##                     'hmGuess':hmGuess,'Dguess':Dguess,'theDiff':theDiff}
    ## if result['X']==0 and theDiff < 1.e-5:
    ##     print "in findDiff, sat: "
    ##     print output % {'r':r,'T':T,'Dtarget':Dtarget,'p':p,'gz':gz,
    ##                     'hmGuess':hmGuess,'Dguess':Dguess,'theDiff':theDiff}
    return theDiff


def dens_line(Dline,r,p,gz):
    """find rt given density
    """
    hlow=300.e3
    hhigh=320.e3
    #print "-"*50,"thermlib"
    #print "start"
    hPlot=optimize.zeros.brenth(findDiffD, hlow, hhigh, (Dline,r,p,gz),xtol=1.e-5);
    ## print "stop"
    ## print "hPlot: ",hPlot
    return hPlot



def findDiffDthetal(thetalGuess,Dtarget,r,p):
    #parcel is just saturated
    try:
        result=t_uos_thetal(thetalGuess,r,p)
    except Exception as x:
        print("class of x =", x.__class__)
        print("something bad in findDiffDthetal II")
        print("here")
        try:
            (type_, value_, traceback_) = sys.exc_info()
        except Exception as x:
            print("caught: ",x.__class__)
        print("here II")
        print("type_ =", type_)
        print("value_ =", value_)
        print("traceback_ =", traceback_)
        for key in dir(traceback_):
            print("traceback_.%s =" % key, eval("traceback_.%s" % key))
    T=result["T"]
    rl=result["RL"]
    rv=result["RV"]
    #print "calc_dens: here is h,rt,T,rl,qv: ",h,rt,T,rl,qv
    #4.3.5 p. 113
    Dguess=p*1.e3*(1. + rl + rv)/(tc.RD*T*(1. + (rv/tc.EPS)))
    theDiff=Dtarget - Dguess
    return theDiff

def dens_line_thetal(Dline,r,p):
    """find thetal given density
    """


    Tlow=230.
    Thigh=330.
    #print "-"*50,"thermlib dens_line_thetal"
    #print "start"
    try:
        thetalPlot=optimize.zeros.brenth(findDiffDthetal, Tlow, Thigh, (Dline,r,p),xtol=1.e-5);
    except:
        print("in thermlib1: caught",(Dline,r,p))
        raise ValueError
    return thetalPlot



def mix(hp,rtp,hi,rti,p,gz,num,smin,smax):
    """in: hp =liquid water static energy for a subcloud air parcel (J/kg)
           rtp=mixing ratio for a subcloud air parcel (kg/kg)
           hi =liquid water static energy for an environmental air parcel (J/kg)
           rti=mixing ratio for an environmental air parcel (kg/kg)
           p=pressure at mixing level (mb)
           gz=geopotential height gz at mixing level(m2/s2)
           num=number of mixtures to be generated
           smin=minimum fraction of environmental air permitted in mixtures
           smax=maximum fraction of environmental air permitted in  mixtures
       out:library containing mixtures properties
           sm=fraction of environmental air in the mixture
           rm=specific total water content in the mixture
           hm=liquid water static energy in the mixture
           tm=temperature (K),
           trom=density temperature (K)
           thmro=virtual potential temperature including liquid water effect (K)
           rlm=specific liquid water content (kg/kg)
           soru=identification of whether or not the mixture is saturated
    """
    Z=np.zeros;
    hm = Z([num],'d'); rm= Z([num],'d');  tm  =Z([num],'d');
    rlm= Z([num],'d'); trom=Z([num],'d'); soru=Z([num],'i');
    throm=Z([num],'d');dens=Z([num],'d')
    mixture={};
    sm=np.linspace(smin, smax, num);
    for j in range(num):
          hm[j]=sm[j]*hi +hp *(1-sm[j]);
          rm[j]=sm[j]*rti+rtp*(1-sm[j]);
          result=t_uos(hm[j],rm[j],p,gz);
          tm  [j] =result["T"];
          rlm [j] =result["RL"];
          soru[j] =result["X"]   
          trom[j] =tro_trtp(tm[j],rm[j],p);
          throm[j]=thv_tvp(trom[j],p); dens[j]=calc_dens(hm[j],rm[j],p,gz)
    mixture["SM"]   =sm;
    mixture["RM"]   =rm;
    mixture["HM"]   =hm;
    mixture["TM"]   =tm;
    mixture["RLM"]  =rlm;
    mixture["SORU"] =soru;  
    mixture["TROM"] =trom;
    mixture["THROM"]=throm;
    mixture["DENS"]=dens

    return mixture

def mix_thetal(thetalp,rtp,thetali,rti,p,num,smin,smax):
    """in: hp =liquid water static energy for a subcloud air parcel (J/kg)
           rtp=mixing ratio for a subcloud air parcel (kg/kg)
           hi =liquid water static energy for an environmental air parcel (J/kg)
           rti=mixing ratio for an environmental air parcel (kg/kg)
           p=pressure at mixing level (mb)
           gz=geopotential height gz at mixing level(m2/s2)
           num=number of mixtures to be generated
           smin=minimum fraction of environmental air permitted in mixtures
           smax=maximum fraction of environmental air permitted in  mixtures
       out:library containing mixtures properties
           sm=fraction of environmental air in the mixture
           rm=specific total water content in the mixture
           hm=liquid water static energy in the mixture
           tm=temperature (K),
           trom=density temperature (K)
           thmro=virtual potential temperature including liquid water effect (K)
           rlm=specific liquid water content (kg/kg)
           soru=identification of whether or not the mixture is saturated
    """
    Z=np.zeros;
    thetalm = Z([num],'d'); rm= Z([num],'d');  tm  =Z([num],'d');
    rlm= Z([num],'d'); trom=Z([num],'d'); soru=Z([num],'i');
    throm=Z([num],'d');dens=Z([num],'d')
    mixture={};
    sm=np.linspace(smin, smax, num);
    for j in range(num):
          thetalm[j]=sm[j]*thetali +thetalp*(1-sm[j]);
          rm[j]=sm[j]*rti+rtp*(1-sm[j]);
          result=t_uos_thetal(thetalm[j],rm[j],p);
          tm  [j] =result["T"];
          rlm [j] =result["RL"];
          soru[j] =result["X"]   
          trom[j] =tro_trtp(tm[j],rm[j],p);
          throm[j]=thv_tvp(trom[j],p); dens[j]=calc_dens_thetal(thetalm[j],rm[j],p)
    mixture["SM"]   =sm;
    mixture["RM"]   =rm;
    mixture["thetal"]   =thetalm;
    mixture["TM"]   =tm;
    mixture["RLM"]  =rlm;
    mixture["SORU"] =soru;  
    mixture["TROM"] =trom;
    mixture["THROM"]=throm;
    mixture["DENS"]=dens

    return mixture


def findDiffT(Tguess,htarget,rtTarget,press,gzval):
    try:
        rsat=rs_tp(Tguess,press)
    except:
        print("class of x =", x.__class__)
        print("something bad in findDiffT")
    CPn=(tc.CPD + rtTarget*tc.CPV)
    gzvalp=(1+rtTarget)*gzval
    rl=rtTarget - rsat
    if rl < 0:
        rl=0.
    if rtTarget > rsat:
        LV=L_t(Tguess)
        hguess=CPn*Tguess - LV*rl + gzvalp
        sat=True
    else:
        hguess=CPn*Tguess + gzvalp
        sat=False
    #print "in findDiffT: ",Tguess,htarget,rtTarget,p,gz,rs,hguess,sat
    theDiff=htarget - hguess
    return theDiff


if __name__=='__main__':
    from . import thermconst as tc
    from .readsound import readsound
    import matplotlib.pyplot as plt
    from .thermlib import LCL

    
    soundDict=readsound("sound.dat")
    gz=soundDict['heightM']*tc.g
    p=soundDict['pkPa']
    Temp=280.
    press=p[50]
    gzval=gz[50]
    rsat=rs_tp(Temp,press)
    rt=11.e-3
    rl=rt - rsat
    if rl < 0:
        rl=0.
    L=L_t(Temp)
    CPn=(tc.CPD + rt*tc.CPV)
    ht=CPn*Temp - L*rl + gzval*(1. + rt)
    Tlow=230.
    Thigh=340.
    theAnswer=optimize.zeros.brenth(findDiffT, Tlow, Thigh, (ht,rt, press,gzval));
    print("temp: ",theAnswer)
    print(t_uos(ht,rt,press,gzval))
    
def thetaes(p, T, ws):
    """in: p=pressure in kPa, T=temperature in Kelvin,
           ws=saturation mixing ratio kg/kg
       out:thetaes=saturated equivalent potential temperature in Kelvin
    """
    thetaes=T*(100./p)**(RD/CPD)*np.exp((LV0*ws)/(CPD*T))
    return thetaes


def find_rs(T, p):
    """in: T=temperature in Kelvin, p=pressure in kPa
       out:ws=saturation water vapor mixing ratio in (kg/kg)
    """
    ws = 0.622*esat(T)/(p-esat(T))
    return ws

def esat(T):
    """in: T=temperature in Kelvin
       out:esat=saturation vapour pressure over water in kPa
    """
    esat = 0.1*10.**(23.832241-5.02808*np.log10(T)-1.3816e-7*10.**(11.344-0.0303998*T)+\
                   8.1328e-3*10.**(3.49149-1302.8844/T)-2949.076/T)
    return esat

def LCL(Tds, Ts, ps):
    """in: Tds=surface dewpoint temperature in Kelvin,
           Ts =surface temperature in Kelvin,
           ps =surface pressure in kPa
       out:LCL=pressure at the LCL in kPa
    """
    wv = find_rs(Tds, ps)
    ao = find_theta(Ts, ps)
    phold = ps
    x, count = 1, 0	
    while abs(x) > 0.01:
        x = .02*(tmr(wv, phold) - tda(ao, phold))
        phold = phold*2.**x
        count = count + 1
        if count == 10:
            break
    LCL=phold
    return LCL

def tmr(ws, p):
    """in: ws=saturation mixing ratio (kg/kg), p=pressure in kPa
       out:tmr=temperature in Kelvin
    """
    
# inverts the CC equation for temperature. First compute saturation
# vapour pressure and then use an approximation to the inverse saturation
# vapour pressure function to compute the temperature.
    ws=1000.*ws  #convert to g/kg
    p=10.*p  #convert to hPa
    x = np.log10(ws*p/(622.+ws))
    tmr = 273.16 + 10.**(.0498646455*x+2.4082965)-280.23475+38.9114*((10.**(.0915*x)-1.2035)**2)
    return tmr

def tda(theTheta, p):
    """in: theTheta=potential temperature in Kelvin, p=pressure in kPa
       out:tda=temperature in Kelvin
    """
    tda =  theTheta*(p/100.)**0.288
    return tda

def thetae(tds, ts, ps):
    """in: Tds=surface dewpoint temperature in Kelvin
           Ts =surface temperature in Kelvin
           ps =surface pressure in kPa
       out:thetae=equivalent potential temperature in Kelvin
    """
    plcl = LCL(tds, ts, ps)
    thetaLcl = find_theta(ts, ps)
    Tlcl = tda(thetaLcl, plcl)
    thetae = find_theta(ts, ps)*np.exp(2.6518986*find_rs(tds, ps)*1000./Tlcl)
    return thetae

def testit():
    "run a default test on the module"
    #from PS 2, problem 2
    p=100.
    T=288.15
    Td=277.15
    print("simple test run: \n")
    print("Temperature (K): ", T)
    print("Dewpoint (K): ",Td)
    print("wv (kg/kg): ", "%5.4f" % find_rs(Td, p))
    print("pressure (kPa): ",p)
    print("pressure LCL (kPa): ","%5.2f" % (LCL(Td,T,p),))
    print("Thetae (K): ","%5.2f" % (thetae(Td,T,p),))
    answer="""
      Temperature (K):  288.15
      Dewpoint (K):  277.15
      wv (kg/kg):  0.0051
      pressure (Pa):  100.0
      pressure LCL (Pa):  84.86
      Thetae (K):  302.67
    """
    print("\nanswer expected: \n",answer)

if __name__== '__main__':
    testit()
    
    
    
