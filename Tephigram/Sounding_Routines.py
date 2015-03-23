from numpy import load, array, radians, sin, cos, linspace, mean, log, isnan, nan, nanmin, nanmax, nanmean, abs, zeros, exp, where,\
                  concatenate, diff
import numpy as np
from numpy.ma import masked_array
from scipy.interpolate import interp1d

# <codecell>

# Exception handling, with line number and stuff
#import linecache
#import sys

import pdb

import imp
imp.load_source('GenMeteoFuncs', '/nfs/see-fs-01_users/eepdw/python_scripts/modules/GeneralMeteoFunctions.py')
from GenMeteoFuncs import *

imp.load_source('MonashFuncs', '/nfs/see-fs-01_users/eepdw/python_scripts/modules/NotYetUsedOrTestedFunctions.py')
from MonashFuncs import *

Rs_da=287.05          # Specific gas const for dry air, J kg^{-1} K^{-1}
Rs_v=461.51           # Specific gas const for water vapour, J kg^{-1} K^{-1}
Cp_da=1004.6          # Specific heat at constant pressure for dry air
Cv_da=719.            # Specific heat at constant volume for dry air
Cp_v=1870.            # Specific heat at constant pressure for water vapour
Cv_v=1410.            # Specific heat at constant volume for water vapour
Cp_lw=4218            # Specific heat at constant pressure for liquid water
Epsilon=0.622         # Epsilon=R_s_da/R_s_v; The ratio of the gas constants
degCtoK=273.15        # Temperature offset between K and C (deg C)
rho_w=1000.           # Liquid Water density kg m^{-3}
grav=9.81             # Gravity, m s^{-2}
Lv=2.5e6              # Latent Heat of vaporisation 
boltzmann=5.67e-8     # Stefan-Boltzmann constant
mv=18.0153            # Mean molar mass of water vapor(g/mol)


# A few basic checks and possible conversions

def PressureinhPaCheck(pres):
    try:
	    maxpres=max(pres)
    except TypeError:
	    maxpres=pres

    if maxpres>2000:
	    print "WARNING: P>2000 hPa; did you input a value in Pa?"
    
# If necessary, convert to celsius

def TempInCelsiusCheck(tc):
    
    if all(tc[~isnan(tc)]>100.):
        tc -= 273.15
        
    return tc

# If necessary, convert to kelvin
    
def TempInKelvinCheck(tk):
    
    if all(tk[~isnan(tk)]<100.):
        tk += 273.15
        
    return tk

# Routines to find enivoronmental temperature/dewpoint at pressure level, based on
# linear interpolation of sounding

def TempC500mb(tempc, pres):
    """Temperature in Celsius at 500 mb

    INPUTS: 
    Temperature profile from sounding - tempc (C)
    Pressure profile from sounding pres (hPa)

    OUTPUTS: Temp(C)

    Source: http://weather.uwyo.edu/upperair/indices.html
    Prints a warning if a pressure value below 2000 Pa input, to ensure
    that the units were input correctly.
    """
    PressureinhPaCheck(pres)
    tempc = TempInCelsiusCheck(tempc)
    
    return interp_point(tempc,pres,500.)

def TempC850mb(tempc, pres):
    """Temperature in Celsius at 500 mb

    INPUTS: 
    Temperature profile from sounding - tempc (C)
    Pressure profile from sounding pres (hPa)

    OUTPUTS: Temp(C)

    Source: http://weather.uwyo.edu/upperair/indices.html
    Prints a warning if a pressure value below 2000 Pa input, to ensure
    that the units were input correctly.
    """
    PressureinhPaCheck(pres)
    tempc = TempInCelsiusCheck(tempc)
    
    return interp_point(tempc,pres,850.)

def DewTempC850mb(dew_tempc, pres):
    """Temperature in Celsius at 500 mb

    INPUTS: 
    Temperature profile from sounding - tempc (C)
    Pressure profile from sounding pres (hPa)

    OUTPUTS: Temp(C)

    Source: http://weather.uwyo.edu/upperair/indices.html
    Prints a warning if a pressure value below 2000 Pa input, to ensure
    that the units were input correctly.
    """

    PressureinhPaCheck(pres)
    dew_tempc = TempInCelsiusCheck(dew_tempc)
    
    return interp_point(dew_tempc,pres,850.)

def TempC700mb(tempc, pres):
    """Temperature in Celsius at 500 mb

    INPUTS: 
    Temperature profile from sounding - tempc (C)
    Pressure profile from sounding pres (hPa)

    OUTPUTS: Temp(C)

    Source: http://weather.uwyo.edu/upperair/indices.html
    Prints a warning if a pressure value below 2000 Pa input, to ensure
    that the units were input correctly.
    """
    PressureinhPaCheck(pres)
    tempc = TempInCelsiusCheck(tempc)
    
    return interp_point(tempc,pres,700.)

def DewTempC700mb(dew_tempc, pres):
    """Temperature in Celsius at 500 mb

    INPUTS: 
    Temperature profile from sounding - tempc (C)
    Pressure profile from sounding pres (hPa)

    OUTPUTS: Temp(C)

    Source: http://weather.uwyo.edu/upperair/indices.html
    Prints a warning if a pressure value below 2000 Pa input, to ensure
    that the units were input correctly.
    """
    PressureinhPaCheck(pres)
    dew_tempc = TempInCelsiusCheck(dew_tempc)
    
    return interp_point(dew_tempc,pres,700.)
   
def MeanFirst500m(vr, height, st_height):
    
    """Average of variable for first 500m above surface

    INPUTS: 
    Sounding variable - vr
    Heights of input sounding = height
    Station height - st_height

    OUTPUTS: Variable average for first 500m above surface

    Source: http://weather.uwyo.edu/upperair/indices.html
    Prints a warning if a pressure value below 2000 Pa input, to ensure
    that the units were input correctly.
    """
    # Calculate average for first 500m
    
    fifty_m_above_surface=st_height+50.
    fivehundred_m_above_surface=st_height+500.
        
    y_points = linspace(fifty_m_above_surface, fivehundred_m_above_surface, 200) # Points for height interpolation in first 500m
    
    #print y_points
    #print vr
    #print height
    
    vr_500m = interp_sounding(vr,height,y_points)
    mean_vr = nanmean(vr_500m)
    
    #print fifty_m_above_surface
    #print fivehundred_m_above_surface
    
    return mean_vr

# Routines to lift parcel from various starting characteristics

def TCParcelLiftedFrom850To500(tempc, dew_tempc, pres):
    """Temperature in Celsius at 500 mb of a parcel lifted from 850 mb 

    INPUTS: 
    Temperature profile from sounding - tempc (C)
    Dewpoint temperature profile from sounding - tempc (C)
    Pressure profile from sounding pres (hPa)

    OUTPUTS: Temp(C)

    Source: http://weather.uwyo.edu/upperair/indices.html
    Prints a warning if a pressure value below 2000 Pa input, to ensure
    that the units were input correctly.
    """
    try:
	    maxpres=max(pres)
    except TypeError:
	    maxpres=pres

    if maxpres>2000:
	    print "WARNING: P>2000 hPa; did you input a value in Pa?"
    
    tempc = TempInCelsiusCheck(tempc)
    
    dew_tempc = TempInCelsiusCheck(dew_tempc)
       
    t850c = interp_point(tempc,pres,850.)
    td850c = interp_point(dew_tempc,pres,850.)
    
    #print t850c
    #print td850c
    
    parcel_profile = ParcelAscentDryToLCLThenMoistC(850.,t850c,td850c, pres)
    
    t500c_lift_from_850 = interp_point(parcel_profile,pres,500.)
    
    return t500c_lift_from_850
    
def TCParcelLiftedFromFirst500mTo500(tempc, dew_tempc, pres, heights, st_height):
    """Temperature in Celsius at 500 mb of a parcel lifted from 850 mb 

    INPUTS: 
    Temperature profile from sounding - tempc (C)
    Dewpoint temperature profile from sounding - tempc (C)
    Pressure profile from sounding pres (hPa)

    OUTPUTS: Temp(C)

    Source: http://weather.uwyo.edu/upperair/indices.html
    Prints a warning if a pressure value below 2000 Pa input, to ensure
    that the units were input correctly.
    """
    try:
	    maxpres=max(pres)
    except TypeError:
	    maxpres=pres

    if maxpres>2000:
	    print "WARNING: P>2000 hPa; did you input a value in Pa?"
    
    tempc = TempInCelsiusCheck(tempc)   
    dew_tempc= TempInCelsiusCheck(dew_tempc)
             
    tc_first_500m = MeanFirst500m(tempc, heights, st_height)
    tdc_first_500m = MeanFirst500m(dew_tempc, heights, st_height)
    p_first_500m = MeanFirst500m(pres, heights, st_height)
    
    parcel_profile = ParcelAscentDryToLCLThenMoistC(p_first_500m,tc_first_500m,tdc_first_500m, pres)
    
    t500c_lift_from_first_500m = interp_point(parcel_profile,pres,500.)
    
    #print t500c_lift_from_first_500m
    
    return t500c_lift_from_first_500m

def LiftDry(startp,startt,startdp,y_points):
    """Lift a parcel to discover certain properties.
    INPUTS:
    startp:  Pressure (hPa)
    startt:  Temperature (C)
    startdp: Dew Point Temperature (C)
    """
    #if startt<startdp:
     #   pdb.set_trace()

    assert startt>startdp,"Not a valid parcel. Check Td<Tc %s %s" % (startt, startdp)
    #Pres=linspace(startp,100,100)

    if startt>100.:
        startt=startt-273.15
    # Lift the dry parcel
    T_dry=(startt+273.15)*(y_points/startp)**(Rs_da/Cp_da)-273.15 
    
    return T_dry

# Routines from https://github.com/tchubb/SkewT/blob/master/build/lib.linux-x86_64-2.7/skewt/thermodynamics.py

def LiftWet(startt,pres):
    #--------------------------------------------------------------------
    # Lift a parcel moist adiabatically from startp to endp.
    # Init temp is startt in C, pressure levels are in hPa    
    #--------------------------------------------------------------------
    if startt>100.:
        startt=startt-273.15
        
    if pres[0]<pres[-1]:
        pres = pres[::-1]
    
    temp=startt
    t_out=zeros(pres.shape);t_out[0]=startt
    for ii in range(pres.shape[0]-1):
	    delp=pres[ii]-pres[ii+1]
 	    temp=temp-100*delp*GammaW(temp+273.15,(pres[ii]-delp/2)*100)
	    t_out[ii+1]=temp

    return t_out[::-1]

# Several indices
# http://weather.uwyo.edu/upperair/indices.html

def KIndex(tempc, dew_tempc, pres):
    T850c = TempC850mb(tempc, pres)
    T500c = TempC500mb(tempc, pres)
    TD850c = DewTempC850mb(dew_tempc,pres)
    T700c = TempC700mb(tempc,pres)
    TD700c = DewTempC700mb(dew_tempc,pres)
    
    return (T850c - T500c) + TD850c - (T700c - TD700c)

def CrossTotalsIndex(tempc,dew_tempc, pres):
    TD850c = DewTempC850mb(dew_tempc,pres)
    T500c = TempC500mb(tempc, pres)
    
    return TD850c - T500c

def VerticalTotalsIndex(tempc,pres):
    T850c = TempC850mb(tempc,pres)
    T500c = TempC500mb(tempc, pres)
    
    return T850c - T500c

def TotalTotalsIndex(tempc, dew_tempc, pres):
    T850c = TempC850mb(tempc,pres)
    T500c = TempC500mb(tempc, pres)  
    TD850c = DewTempC850mb(tempc,pres)
    
    return (T850c - T500c) + (TD850c - T500c)

def LiftedIndex(tempc, dew_tempc, pres, heights, st_height):
    """LIFT	= T500 - Tparcel
		T500	= temperature in Celsius of the environment at 500 mb
		Tparcel	= 500 mb temperature in Celsius of a lifted parcel with 
        the average pressure, temperature, and dewpoint of the layer 500 m 
        above the surface 

    INPUTS: 
    500mb temp of parcel lifted from 850mb (C)
    500mb temp (C)
    pref: 

    OUTPUTS: Temp(C)

    Source: http://glossary.ametsoc.org/wiki/Stability_index
            http://weather.uwyo.edu/upperair/indices.html
    Prints a warning if a pressure value below 2000 Pa input, to ensure
    that the units were input correctly.
    """
    T500c = TempC500mb(tempc, pres)
    
    
    t500c_lift_from_first_500m = TCParcelLiftedFromFirst500mTo500(tempc, dew_tempc, pres, heights, st_height)
    
    #print t500c_lift_from_first_500m
    #print T500c
        
    lift = T500c - t500c_lift_from_first_500m
    
    return lift

def ShowalterIndex(tempc, dew_tempc, pres):
    """SHOW	= T500 - Tparcel
		T500	= Temperature in Celsius at 500 mb
		Tparcel	= Temperature in Celsius at 500 mb of a parcel lifted from 850 mb 

    INPUTS: 
    500mb temp of parcel lifted from 850mb (C)
    500mb temp (C)
    pref: 

    OUTPUTS: Temp(C)

    Source: http://glossary.ametsoc.org/wiki/Stability_index
    Prints a warning if a pressure value below 2000 Pa input, to ensure
    that the units were input correctly.
    """
    t500c = TempC500mb(tempc,pres)
    
    t500c_lift_from_850 = TCParcelLiftedFrom850To500(tempc, dew_tempc, pres)
        
    show = t500c - t500c_lift_from_850
    
    return show

#Interpolate Sounding

#y_points=linspace(5000, 100000, 200) # Points for pressure interpolation
    
def interp_sounding(variable, pressures_s,y_points):
    
    #print y_points
    #print variable

    #pdb.set_trace()

    #nan_mask = masked_array(array(variable, dtype=float), isnan(array(variable, dtype=float)))
    #nan_mask_p = masked_array(array(pressures_s, dtype=float), isnan(array(variable, dtype=float)))
  try:
     y_interp = np.empty(y_points.shape)
     y_interp.fill(np.nan)

     if type(pressures_s)!=float and type(variable)!=float:
 
        if len(np.array(pressures_s))>1 and len(np.array(variable))>1:

            variable = [x for (y,x) in sorted(zip(pressures_s, variable), key=lambda pair: pair[0])]
            pressures_s = [y for (y,x) in sorted(zip(pressures_s, variable), key=lambda pair: pair[0])]
  
            interp = interp1d(pressures_s, variable, bounds_error=False, fill_value=nan)
            y_interp = interp(y_points)

  except Exception,e:
        pdb.set_trace()
        print e
        print variable
        print pressures_s

  return y_interp

def interp_point(variable, pressures_s,y_points):
    
    #print y_points
    #print variable

    #pdb.set_trace()

    #nan_mask = masked_array(array(variable, dtype=float), isnan(array(variable, dtype=float)))
    #nan_mask_p = masked_array(array(pressures_s, dtype=float), isnan(array(variable, dtype=float)))
  try:

     if type(pressures_s)!=float and type(variable)!=float:
 
        if len(np.array(pressures_s))>1 and len(np.array(variable))>1:

            variable = [x for (y,x) in sorted(zip(pressures_s, variable), key=lambda pair: pair[0])]
            pressures_s = [y for (y,x) in sorted(zip(pressures_s, variable), key=lambda pair: pair[0])]
  
            interp = interp1d(pressures_s, variable, bounds_error=False, fill_value=nan)
            y_interp = interp(y_points)

  except Exception,e:
        pdb.set_trace()
        print e
        print variable
        print pressures_s

  return y_interp
def LiftedCondensationLevelTemp(init_temp_k, dew_init_temp_k): 

    """
    """
    if (init_temp_k<100.):
        init_temp_k = init_temp_k +273.15
    if (dew_init_temp_k<100.):
        dew_init_temp_k = dew_init_temp_k +273.15

    return (1./(1./(dew_init_temp_k-56) + log(init_temp_k/dew_init_temp_k)/800.)) + 56



def LiftedCondensationLevelPres(pres_parcel, temp_parcel_k, dew_temp_parcel_k):

    """
    """

    vap_press = VapourPressure(dew_temp_parcel_k-273.15)
    sat_temp_k = 55+2840/(3.5*log(temp_parcel_k)-log(vap_press/100)-4.805)

    wvmr = WaterVapourMixingRatio(vap_press,pres_parcel)

    return pres_parcel*(sat_temp_k/temp_parcel_k)**((Cp_da+Cp_v*wvmr)/(Rs_da*(1+wvmr/Epsilon)))


def SurfaceMinus(pressures, surface_pressure, pressure_diff):

    '''
    Find closest point to surface_pressure in pressures  minus pressure_diff 
    (Pa)
    if it is not too far away (e.g in a sounding where there are no interpolated values 
    near the surface because of no data in the sounding near the surface)
    '''
          
    # If there is a point not too far away from (surface pressure - 5 hPa)

    if np.abs(pressures-(surface_pressure-pressure_diff)).min()<(700):   
            
         # Closest point to 5 hPa less than surface

         i_idx = np.abs(pressures-(surface_pressure-pressure_diff)).argmin()  

         #if np.isnan(sat_temp_interp[i_idx]):
          #      print sat_temp_interp
           #     st=sat_temp_interp
            #    vp=vap_press_interp
             #   i_d=i_idx
              #  nm=nan_mask
               # #sf=surf_level
                #sp=surface_p

    else:

         i_idx=np.nan

    return i_idx


def PBLFromParcelVPT(pressures, temps_cent, dewp_temps_cent, surface_pressure, pressure_diff_from_surface):

    """
    First method from http://dx.doi.org/10.1029/2009JD013680
    The parcel method [Holzworth, 1964; Seibert et al., 2000]
    """

    #pdb.set_trace()


    pbl_pressure=np.nan

    if (any(temps_cent[~np.isnan(temps_cent)])) & (any(dewp_temps_cent[~np.isnan(dewp_temps_cent)])) & (~np.isnan(surface_pressure)): #If arrays are not totally nan

        nan_mask = np.ma.masked_array(np.array(pressures), np.isnan(np.array(dewp_temps_cent, dtype=float)) | np.isnan(np.array(temps_cent, dtype=float)))

        #pdb.set_trace()

        if (np.sum(nan_mask>50000)>=10) &  (np.abs(nan_mask-(surface_pressure-pressure_diff_from_surface)).min()<(700)):   

            #pdb.set_trace()

            i_idx = SurfaceMinus(nan_mask, surface_pressure, pressure_diff_from_surface)  

            # Calculate virtual potential temperature of profile
   
            environment_vpt = ThetaV(nan_mask, temps_cent, dewp_temps_cent)
   
            surface_parcel_vpt = environment_vpt[i_idx]

            # Mask out surface point in virtual potential temperature profile
            mask = np.zeros( environment_vpt.shape, dtype=bool) # all elements False
            mask[i_idx] = True  # Mask out surface point
            environment_vpt_mask = np.ma.masked_array(environment_vpt, mask)
                               # Find index of point closest to minimum difference
 
                               #   vpt_idx=np.nanargmin(np.abs(parcel_vpt-surface_vpt[i_idx]))  
        
            # Interpolate difference of virtual potential temperature profile from surface vpt and find where equal to zero
    
            interp = interp1d((environment_vpt_mask-surface_parcel_vpt),nan_mask , bounds_error=False, fill_value=np.nan)
     
            pbl_pressure = float(interp(0.))

            #if ~np.isnan(pbl_pressure):
                #pdb.set_trace()

  
    return pbl_pressure



def SatVap(tempc,phase="liquid"):
    """Calculate saturation vapour pressure over liquid water and/or ice.

    INPUTS: 
    tempc: (C)
    phase: ['liquid'],'ice'. If 'liquid', do simple dew point. If 'ice',
    return saturation vapour pressure as follows:

    Tc>=0: es = es_liquid
    Tc <0: es = es_ice

   
    RETURNS: e_sat  (Pa)
    
    SOURCE: http://cires.colorado.edu/~voemel/vp.html (#2:
    CIMO guide (WMO 2008), modified to return values in Pa)
    
    This formulation is chosen because of its appealing simplicity, 
    but it performs very well with respect to the reference forms
    at temperatures above -40 C. At some point I'll implement Goff-Gratch
    (from the same resource).
    """

    over_liquid=6.112*exp(17.67*tempc/(tempc+243.12))*100.
    over_ice=6.112*exp(22.46*tempc/(tempc+272.62))*100.
    # return where(tempc<0,over_ice,over_liquid)

    if phase=="liquid":
	# return 6.112*exp(17.67*tempc/(tempc+243.12))*100.
	return over_liquid
    elif phase=="ice":
	# return 6.112*exp(22.46*tempc/(tempc+272.62))*100.
	return where(tempc<0,over_ice,over_liquid)
    else:
	raise NotImplementedError

def GammaW(tempk,pres,e=None):
    """Function to calculate the moist adiabatic lapse rate (deg C/Pa) based
    on the temperature, pressure, and rh of the environment.

    INPUTS:
    tempk (K)
    pres (Pa)
    RH (%)

    RETURNS:
    GammaW: The moist adiabatic lapse rate (Dec C/Pa)
    """

    tempc=tempk-degCtoK
    es=SatVap(tempc)
    ws=WaterVapourMixingRatio(es,pres)

    if e is None:
	# assume saturated
	e=es

    w=WaterVapourMixingRatio(e,pres)

    tempv=VirtualTempFromMixR(tempk,w)
    latent=Latentc(tempc)

    A=1.0+latent*ws/(Rs_da*tempk)
    B=1.0+Epsilon*latent*latent*ws/(Cp_da*Rs_da*tempk*tempk)
    Rho=pres/(Rs_da*tempv)
    Gamma=(A/B)/(Cp_da*Rho)
    return Gamma

def LFCParcelAscent(parcel_profile, temp_k, pres):

        """

        """

        lfc_idx = nanmax(where((parcel_profile>temp_k) & (pres<(nanmax(pres)-50)))[0])
    
        #lfc_temp = temp_k[lfc_idx]
        lfc_pres = pres[lfc_idx]    
        # If parcel unstable throughout sounding set LFC to LCL
    
        if all(parcel_profile>temp_k):
            #lfc_temp = lcl_temp
            lfc_pres = lcl_pres
            
        return lfc_pres


def EQLVParcelAscent(parcel_profile, temp_k, pres):
    
    temp_k = TempInKelvinCheck(temp_k)
    parcel_profile = TempInKelvinCheck(parcel_profile)
    
    idx_saddles=[]
    for i in (where((parcel_profile>temp_k))[0]):
        if i in (where(parcel_profile<temp_k)[0]+1):
            #print i
            idx_saddles.append(i)
            
    idx_eqlv = nanmin(idx_saddles)
    #print idx_eqlv
    
    # Interpolate around index of highest saddle point to zero difference
    
    y_points_zero=0 # Points from original sounding interpolation
    eqlv_p = interp_point(pres[idx_eqlv-1:idx_eqlv+1],(parcel_profile-temp_k)[idx_eqlv-1:idx_eqlv+1],y_points_zero)
    eqlv_t = interp_point(temp_k[idx_eqlv-1:idx_eqlv+1],(parcel_profile-temp_k)[idx_eqlv-1:idx_eqlv+1],y_points_zero)
        
    return eqlv_p, eqlv_t


def ParcelAscentDryToLCLThenMoistC(init_parcel_pres,init_parcel_temp_c,init_parcel_dew_temp_c, pres):
        dry_parcel_TC = LiftDry(init_parcel_pres,init_parcel_temp_c,init_parcel_dew_temp_c, pres)
        dry_parcel_TK = dry_parcel_TC+273.15
        
        if (init_parcel_temp_c>100.):
            init_parcel_temp_c = init_parcel_temp_c - 273.15
        if (init_parcel_dew_temp_c>100.):
            init_parcel_dew_temp_c = init_parcel_dew_temp_c - 273.15
                 
        lcl_temp = LiftedCondensationLevelTemp(init_parcel_temp_c+273.15, init_parcel_dew_temp_c+273.15) 
            
    
        i_idx = (abs(dry_parcel_TK-lcl_temp)).argmin() 
    
        #  Find initial temp at LCL for moist parcel ascent
    
        y_points_zero=0
    
        moist_t_init = dry_parcel_TC[i_idx]
    
        temp_moist_adi_parcel_above_lcl_c = LiftWet(moist_t_init,pres[0:i_idx])
        #temp_moist_adi_parcel_above_lcl_k = temp_moist_adi_parcel_above_lcl_c + 273.15

        # Combine dry parcel ascent below LCL and moist above

       

        #print temp_moist_adi_parcel_above_lcl_c
  
        parcel_profile=concatenate((temp_moist_adi_parcel_above_lcl_c, dry_parcel_TC[i_idx::]))
        
        return parcel_profile
        

# CAPE and CIN

def CapeCinPBLInput(pres, temp_k, dew_temp_k, height, st_height, pbl_pressure): # PBL_pressure input temporary ?
       
    # Calculate pressure, temperature and dew point temperature averages for first 500m

    
    temp_k = TempInKelvinCheck(temp_k)
    dew_temp_k = TempInKelvinCheck(dew_temp_k)
    
    mean_temp_k = MeanFirst500m(temp_k, height, st_height)
    dew_mean_temp_k = MeanFirst500m(dew_temp_k, height, st_height)
    mean_pres = MeanFirst500m(pres, height, st_height)
    
    #temp_c=TempInKelvinCheck(temp_k)
    #dew_temp_c=TempInKelvinCheck(dew_temp_k)  

    mean_temp_c = mean_temp_k - 273.15
    dew_mean_temp_c = dew_mean_temp_k - 273.15
    
    # Find LCL temp and pressure
    
    lcl_t = LiftedCondensationLevelTemp(mean_temp_k, dew_mean_temp_k)
      
    vap_press = VapourPressure(dew_mean_temp_c)

    wvmr = WaterVapourMixingRatio(vap_press,mean_pres*100)
      
    kappa=PoissonConstant(wvmr)
    
    lcl_p = LiftedCondensationLevelPres(mean_pres, lcl_t, mean_temp_k)
    
    # Calculate dry parcel ascent 
       
    parcel_profile = ParcelAscentDryToLCLThenMoistC(mean_pres,mean_temp_c,dew_mean_temp_c, pres)
    parcel_profile = parcel_profile+273.15
    # Find equilibrium level

    eqlv_p, eqlv_t = EQLVParcelAscent(parcel_profile, temp_k, pres)
    
    # Find LFC
        
    lfc_p = LFCParcelAscent(parcel_profile, temp_k, pres)
    
    # Calculate CAPE and CIN
    
    delta_z=diff(height)    # Find delta height in sounding
    
    # Take all but lowest pressure level (so length matches delta_z)
    
    pp_diff = parcel_profile[1::]
    Tk_diff = temp_k[1::] 
    p_diff=pres[1::]

    #pdb.set_trace()

    sum_ascent = abs(delta_z)*(pp_diff-Tk_diff)/Tk_diff
    
    CAPE = grav*sum(sum_ascent[((pp_diff-Tk_diff)>0) & (p_diff>eqlv_p) & (p_diff<lfc_p)])
    
    #Taking levels above lcl but uwyo specifies top of mixed layer

    #CIN = grav*sum(sum_ascent[((pp_diff-Tk_diff)<0) & (p_diff>lfc_p) & (p_diff<lcl_p)])
    CIN = grav*sum(sum_ascent[((pp_diff-Tk_diff)<0) & (p_diff>lfc_p) & (p_diff<pbl_pressure)])
    
    #print "LCL_Test %s" % lcl_test
    
    #print "Mean Pressure 500m %s" % mean_pres
    #print "Dew Mean Temp 500m C %s" %dew_mean_temp_c
    #print "Mean Temp 500m C %s" %mean_temp_c
    #print "Vapour Pressure %s" %vap_press   
    #print "LCL Pressure %s" % lcl_p
    #print "LFC Pressure %s" % lfc_p
    #print "EQLV Pressure %s" % eqlv_p
    #print "CAPE %s" % CAPE
    #print "CIN %s" % CIN
    
      
    return eqlv_p, parcel_profile,lcl_p, lfc_p, lcl_t, delta_z, CAPE, CIN
# CAPE and CIN

def CapeCin(pres, temp_k, dew_temp_k, height, st_height):
       
    # Calculate pressure, temperature and dew point temperature averages for first 500m

    
    temp_k = TempInKelvinCheck(temp_k)
    dew_temp_k = TempInKelvinCheck(dew_temp_k)
    
    mean_temp_k = MeanFirst500m(temp_k, height, st_height)
    dew_mean_temp_k = MeanFirst500m(dew_temp_k, height, st_height)
    mean_pres = MeanFirst500m(pres, height, st_height)
    
    #temp_c=TempInKelvinCheck(temp_k)
    #dew_temp_c=TempInKelvinCheck(dew_temp_k)  

    mean_temp_c = mean_temp_k - 273.15
    dew_mean_temp_c = dew_mean_temp_k - 273.15
    
    # Find LCL temp and pressure
    
    lcl_t = LiftedCondensationLevelTemp(mean_temp_k, dew_mean_temp_k)
      
    vap_press = VapourPressure(dew_mean_temp_c)

    wvmr = WaterVapourMixingRatio(vap_press,mean_pres*100)
      
    kappa=PoissonConstant(wvmr)
    
    lcl_p = LiftedCondensationLevelPres(mean_pres, lcl_t, mean_temp_k)
    
    # Calculate dry parcel ascent 
       
    parcel_profile = ParcelAscentDryToLCLThenMoistC(mean_pres,mean_temp_c,dew_mean_temp_c, pres)
    parcel_profile = parcel_profile+273.15
    # Find equilibrium level

    eqlv_p, eqlv_t = EQLVParcelAscent(parcel_profile, temp_k, pres)
    
    # Find LFC
        
    lfc_p = LFCParcelAscent(parcel_profile, temp_k, pres)
    
    # Calculate CAPE and CIN
    
    delta_z=diff(height)    # Find delta height in sounding
    
    # Take all but lowest pressure level (so length matches delta_z)
    
    pp_diff = parcel_profile[1::]
    Tk_diff = temp_k[1::] 
    p_diff=pres[1::]

    sum_ascent = abs(delta_z)*(pp_diff-Tk_diff)/Tk_diff
    
    CAPE = grav*sum(sum_ascent[((pp_diff-Tk_diff)>0) & (p_diff>eqlv_p) & (p_diff<lfc_p)])
    
    #Taking levels above lcl but uwyo specifies top of mixed layer

    CIN = grav*sum(sum_ascent[((pp_diff-Tk_diff)<0) & (p_diff>lfc_p) & (p_diff<lcl_p)])
    #CIN = grav*sum(sum_ascent[((pp_diff-Tk_diff)<0) & (p_diff>pbl_pressure) & (p_diff<pbl_pressure)])
    
    #print "LCL_Test %s" % lcl_test
    
    #print "Mean Pressure 500m %s" % mean_pres
    #print "Dew Mean Temp 500m C %s" %dew_mean_temp_c
    #print "Mean Temp 500m C %s" %mean_temp_c
    #print "Vapour Pressure %s" %vap_press   
    #print "LCL Pressure %s" % lcl_p
    #print "LFC Pressure %s" % lfc_p
    #print "EQLV Pressure %s" % eqlv_p
    #print "CAPE %s" % CAPE
    #print "CIN %s" % CIN
    
      
    return eqlv_p, parcel_profile,lcl_p, lfc_p, lcl_t, delta_z, CAPE, CIN

