import numpy as np

Rs_da=287.05          # Specific gas const for dry air, J kg^{-1} K^{-1}
Rs_v=461.51           # Specific gas const for water vapour, J kg^{-1} K^{-1}
Cp_da=1004.6          # Specific heat at constant pressure for dry air
Cv_da=719.            # Specific heat at constant volume for dry air
Cp_v=1870.            # Specific heat at constant pressure for water vapour
Cv_v=1410.            # Specific heat at constant volume for water vapour
Cp_lw=4218            # Specific heat at constant pressure for liquid water
Epsilon=0.622         # Epsilon=R_s_da/R_s_v The ratio of the gas constants
degCtoK=273.15        # Temperature offset between K and C (deg C)
rho_w=1000.           # Liquid Water density kg m^{-3}
grav=9.81             # Gravity, m s^{-2}
Lv=2.5e6              # Latent Heat of vaporisation 
boltzmann=5.67e-8     # Stefan-Boltzmann constant
mv=18.0153            # Mean molar mass of water vapor(g/mol)

def UVWinds(wind_direction, wind_speed):

    """
    """

    wind_rad = np.radians(wind_direction)
    u_wind=-((wind_speed)*np.sin(wind_rad))
    v_wind=-((wind_speed)*np.cos(wind_rad))

    return u_wind,v_wind   

def VapourPressure(dewp_temps_cent):

    """
    Takes dewpoint temperature in centigrade (use temperature for saturated vapour pressure
    and returns saturation vapour pressure
    """

    return 611.2*np.exp(17.67*(dewp_temps_cent/(dewp_temps_cent+243.5)))

def RelativeHumidity(temps_cent, dewp_temps_cent):

    """
    Take temperature and dewpoint in C and returns relative humidity
    """

    #assert 

    sat_vap_pres = VapourPressure(temps_cent)
         
    vap_press = VapourPressure(dewp_temps_cent)

   

    return 100.*vap_press/sat_vap_pres

def WaterVapourMixingRatio(pressures, dewp_temps_cent):

    """
    Takes pressure (Pa) and and dewpoint temperature (C) 
    and returns WVMR (kg/kg) 
    """

    vap_press = VapourPressure(dewp_temps_cent)

    return 0.622*vap_press/(pressures - vap_press) 

def PoissonConstant(wvmr):
    
        """http://glossary.ametsoc.org/wiki/Poisson_constant
        May need to tweak low limit for dry air (=0.2854)"""
     
        PoC = np.where(wvmr>0., 0.2854*(1-0.24*wvmr), 0.2854)
        PoC = np.where(~np.isnan(wvmr), PoC, np.nan)
        return PoC



def SpecificHumidity(pressures, dewp_temps_cent):

    """

    """

    wvmr = WaterVapourMixingRatio(pressures, dewp_temps_cent)

    return wvmr/(1+wvmr)

def Theta(pressures, temps_cent, dewp_temps_cent):

    """

    """

    temp_k = temps_cent + 273.15

    wvmr = WaterVapourMixingRatio(pressures, dewp_temps_cent)
        
    kappa = PoissonConstant(wvmr)

    pressures_hpa = np.array(pressures)/100.

    return temp_k*((1000./pressures_hpa)**(kappa))

# def Theta(temp_k,pres,pref=100000.):
#     """Potential Temperature

#     INPUTS: 
#     temp_k (K)
#     pres (Pa)
#     pref: 

#     OUTPUTS: Theta (K)

#     Source: Wikipedia
#     Prints a warning if a pressure value below 2000 Pa input, to ensure
#     that the units were input correctly.
#     """

#     return temp_k*(pref/pres)**(Rs_da/Cp_da)

def ThetaE(pressures, temps_cent, dewp_temps_cent):

    """

    """

    temp_k = temps_cent + 273.15

    pressures_hpa = np.array(pressures)/100.

    wvmr = WaterVapourMixingRatio(pressures, dewp_temps_cent)

    vap_press =  VapourPressure(dewp_temps_cent)

    sat_temp = 55.+2840./(3.5*np.log(temps_cent)-np.log(vap_press/100.)-4.805)

    return temp_k*((1000./pressures_hpa)**(0.2854*(1.-0.28*wvmr)))*np.exp(((3376./sat_temp)-2.54)*wvmr*(1.+0.81*wvmr))

def ThetaV(pressures, temps_cent, dewp_temps_cent):
    """Virtual Potential Temperature
    
    INPUTS
    temp_k (K)
    pres (Pa)
    e: Water vapour pressure (Pa) (Optional)
    """ 

    mixr = WaterVapourMixingRatio(pressures, dewp_temps_cent)
    theta = Theta(pressures, temps_cent, dewp_temps_cent)

    return theta*(1+mixr/Epsilon)/(1+mixr)

def SaturationTemperature(temps_cent, dewp_temps_cent):

    """

    """

    temp_k = temps_cent + 273.15

    vap_press = VapourPressure(dewp_temps_cent)

    return 55.+2840./(3.5*np.log(temps_cent)-np.log(vap_press/100.)-4.805) 
