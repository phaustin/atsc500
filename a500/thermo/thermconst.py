CPD=1004.5     #heat capacity of dry air (J/kg/K)
CPV=1870.0     #heat capacity of water vapor (J/kg/K)
CL=4190.0      #heat capacity of liquid water (J/kg/K)
RV=461.50      #gas constant for water vapor (J/kg/K)
RD=287.04      #gas constant for dry air (J/kg/K)
LV0=2.501E6    #Latent heat of vaporization at T=273.15(J/kg)
g=9.80616      #gravity acceleration (m/s2)
T0=273.16      #water triple point temperature Kelvin (K)
e0=0.6112      #water triple point pressure (kPa)
P0=100         #reference pressure for potential temperature(kPa)
CPVMCL=CL-CPV  #derived constant
EPS=RD/RV      #derived constant
RDOCP=RD/CPD   #derived constant
EPSI=1./EPS - 1.    #derived constant for Tv
