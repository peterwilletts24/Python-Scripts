  #Raw data read - Climatology
    # /nfs/a90/eepdw/Data/Observations/Radiosonde_downloaded_from_NOAA_GUAN/igra-stations.txt
    # /nfs/a90/eepdw/Data/Observations/Radiosonde_downloaded_from_NOAA_GUAN/readme.txt
    
    # Station numbers for Bombay - 43003; New Delhi - 42182, Lucknow - 42369, Cochin, Bangalore, Madras, Minicoy, Nagpur Sonegaon)
import datetime
import scipy.interpolate
import re
import numpy as np
    
import os 
    
import imp

imp.load_source('SoundingRoutinesTemp', '/nfs/see-fs-01_users/eepdw/python_scripts/Tephigram/Sounding_Routines_Mess.py')
from SoundingRoutinesTemp import *

imp.load_source('GenMeteoFuncs', '/nfs/see-fs-01_users/eepdw/python_scripts/modules/GeneralMeteoFunctions.py')
from GenMeteoFuncs import *
imp.load_source('SoundingRoutines', '/nfs/see-fs-01_users/eepdw/python_scripts/Tephigram/Sounding_Routines.py')
from SoundingRoutines import *
imp.load_source('GeogFuncs', '/nfs/see-fs-01_users/eepdw/python_scripts/modules/GeogFunctions.py')
from GeogFuncs import *
station_list_cs=[43003, 43014, 42867, 43371, 43353, 43285,  43192, 43150, 42339, 40990, 40948]

y_points=np.linspace(5000, 108000, 200) # Points for pressure interpolation

#match_header = re.compile(r'(#.....20..0[56789|#.....19..0[56789])')
match_header = re.compile(r'(#.....20|#.....19)')
#match_header_20th = re.compile(r'#.....19..0[56789]')
match_bad = re.compile(r'9999')


variable_list={'pressures': 0, 'temps':1, 'dewpoints':2, 'winddirs':3, 'windspeeds':4, 'pot_temp':5, 
               'sat_vap_pres':6, 'vap_press':7, 'rel_hum':8, 'wvmr':9, 'sp_hum':10, 'sat_temp':11, 
               'theta_e':12, 'theta_e_sat':13, 'theta_e_minus_theta_e_sat':14, 'geop_height':15, 'u_wind':16, 'v_wind':17}
variable_list_line={'lcl_temp': 0, 'lcl_vpt':1, 'pbl_pressure':2, 'surface_pressure':3, 'T_eq_0':4, 
                    'CAPE':5, 'CIN':6, 'lclp':7, 'lclt':8, 'lfcp':9, 'equil_level':10}
    
test_count=0
    
    
for stat in station_list_cs:
        
      i = '/nfs/a90/eepdw/Data/Observations/Radiosonde_downloaded_from_NOAA_GUAN/%s.dat' % stat

      st,la,lo, st_height = StationInfoSearch(stat)
        
      lev_count=False
        
      single_values=[]
      pressures=[]
      temps=[]
      dewpoints=[]
      winddirs=[]
      windspeeds=[]
      geop_height=[]
        
      pot_temp=[]
      sat_vap_pres=[]
      vap_press=[]
      rel_hum=[]
      wvmr=[]
      sp_hum=[]
      sat_temp=[]
      theta_e=[]
      theta_e_sat=[]
      theta_es_min_theta_e=[]
      temp_lcl=[]
      virtual_pot_temp=[]
      
      #lcl=[]

      pbl_pressure_l=[]
      
      surface_p_l=[]
        
      t_zero_line=[]
        
      CAPE=[]
      CIN=[]
      lclp=[]
      lclt=[] 
      lfcp=[]
      equil_level=[]
 
      #temp_lcl_s=[]
      #virtual_pot_temp_s=[]
      #pressures_for_plotting = np.array((dates_for_plotting, pressures, temps, dewpoints, winddirs, windspeeds))
        
      dates_for_plotting_single=[]
      dates_for_plotting_single_single=[]
        
      levs=-1
      numlevs=0
       
      print i

      f = open(i,'r')
      #print f
      for line in f:
            #line = line.strip().replace('-99999',' %s ' %float('NaN'))
            line = line.strip()
            #print line
            #columns=line.split()
                
            # ###### Extract independent headers #########
    
            #print uwind
            #if match_header_20th.search(line) is not None or match_header_21st.search(line) is not None:
            if match_header.search(line) is not None:
                #print line
        
                year=int(line[6:10])
                month=int(line[10:12])
                day=int(line[12:14])
                hour=int(line[14:16])
        
                date = datetime.datetime(year,month,day,hour)
        
                numlevs=int(line[20:24])
        
                #print numlevs
        
                lev_count=True
                levs=0
        
                pressures_s=[]
                temps_s=[]
                dewpoints_s=[]
                winddirs_s=[]
                windspeeds_s=[]
                geop_height_s=[]

            if (lev_count is True) & (match_header.search(line) is None): 
        
                 # Line length varies depending on length of first variable (pressure)
                 #   li=len(line)
                 # ld=142-len(line)
                 #print line[0:2]

                 #pdb.set_trace()

                if float(line[0:2])==21:
                      #pdb.set_trace()

                      try:
                           surface_p=float(line[2:8].replace('-9999','%s' %float('NaN')).replace('-8888','%s' %float('NaN')))
                           surface_p_bool=True
                           #print surface_p
                      except Exception,e:
                           print e
                           print line
                try:
                    pressures_s.append(float(line[2:8].replace('-9999','%s' %float('NaN')).replace('-8888','%s' %float('NaN'))))     
                    temps_s.append(float(line[15:20].replace('-9999','%s' %float('NaN')).replace('-8888','%s' %float('NaN')))) 
                    dewpoints_s.append(float(line[21:26].replace('-9999','%s' %float('NaN')).replace('-8888','%s' %float('NaN')))) 
                    winddirs_s.append(float(line[26:31].replace('-9999','%s' %float('NaN')).replace('-8888','%s' %float('NaN')))) 
                    windspeeds_s.append(float(line[31:36].replace('-9999','%s' %float('NaN')).replace('-8888','%s' %float('NaN'))))
                    geop_height_s.append(float(line[9:14].replace('-9999','%s' %float('NaN')).replace('-8888','%s' %float('NaN'))))
         
       
                    levs+=1
                except Exception,e:
                    print e
                    print line


            if (match_header.search(line) is None) and levs >= numlevs and pressures_s and list(temps_s) and list(dewpoints_s) and winddirs_s and windspeeds_s:
                    
                 lev_count=False
        
                 if ~np.all(np.isnan(pressures_s)) and ~np.all(np.isnan(temps_s)) and ~np.all(np.isnan(dewpoints_s))\
                    and ~np.all(np.isnan(winddirs_s)) and ~np.all(np.isnan(windspeeds_s))\
                    and len(pressures_s)>1 and len(temps_s)>1 and len(dewpoints_s)>1:
        
                      # Calculate sounding variables

                      #print pressures_s
                      #print temps_s
        
                      #temps_s_mask = np.ma.masked_array(np.array(temps_s, dtype=float), np.isnan(np.array(temps_s, dtype=float)))
                      #pressures_s_mask =  np.ma.masked_array(np.array(pressures_s, dtype=float), np.isnan(np.array(pressures_s, dtype=float)))   
                      #dewpoints_s_mask = np.ma.masked_array(np.array(dewpoints_s, dtype=float), np.isnan(np.array(dewpoints_s, dtype=float)))
                      geop_height_s = np.ma.masked_array(np.array(geop_height_s, dtype=float), np.isnan(np.array(geop_height_s, dtype=float)))
           
                      temps_scent = np.array(temps_s, dtype=float)/10.
                      temps_s = np.array(temps_s, dtype=float)/10. + 273.15
                      pressures_s_hPa = np.array(pressures_s, dtype=float)/100. # Pressures in raw files as mb*100 - which is Pa - convert to hPa for calcs
                      dewpoints_s = np.array(dewpoints_s, dtype=float)/10.
            
                      dewpoint_temp_cent = temps_scent-dewpoints_s
        
                      sat_vap_pres_s = VapourPressure(temps_scent)
         
                      vap_press_s = VapourPressure(dewpoint_temp_cent)
            
                      rel_hum_s = RelativeHumidity(temps_scent, dewpoint_temp_cent)
    
                      wvmr_s = WaterVapourMixingRatio(pressures_s, dewpoint_temp_cent)
                      wvmr_sat_s = WaterVapourMixingRatio(pressures_s, temps_scent)

                      pot_temp_s = Theta(pressures_s, temps_scent, dewpoint_temp_cent)
            
                      sp_hum_s = SpecificHumidity(pressures_s, dewpoint_temp_cent)
            
                      sat_temp_s = SaturationTemperature(temps_scent, dewpoint_temp_cent)
            
                      sat_temp_sat_s = SaturationTemperature(temps_scent, dewpoint_temp_cent)
        
                      theta_e_s = ThetaE(pressures_s, temps_scent, dewpoint_temp_cent)           
      
                      theta_e_sat_s = ThetaE(pressures_s, temps_scent, temps_scent) 
         
    
                      # Interpolate sounding variables
                      #try:
                      temp_interp = interp_sounding(temps_s, pressures_s,y_points)
                
                      #print pressures_s
                      dewpoints_interp = interp_sounding(dewpoints_s, pressures_s,y_points)
                      #print dewpoints_interp
                      winddirs_interp= interp_sounding(np.array(winddirs_s, dtype=float), pressures_s, y_points)
                      windspeeds_interp= interp_sounding(np.array(windspeeds_s, dtype=float)/10, pressures_s, y_points)
                      #geop_height_interp=interp_sounding(np.array(geop_height_s, dtype=float), pressures_s, y_points)
                      geop_height_interp=interp_sounding(geop_height_s.data[~geop_height_s.mask], np.array(pressures_s)[~geop_height_s.mask], y_points)
                      
                      pot_temp_interp = interp_sounding(pot_temp_s, pressures_s,y_points)
                      #print pot_temp_interp
                      sat_vap_pres_interp = interp_sounding(sat_vap_pres_s, pressures_s,y_points)
                      vap_press_interp = interp_sounding(vap_press_s, pressures_s,y_points)
                      rel_hum_interp = interp_sounding(rel_hum_s, pressures_s,y_points)
                      #print rel_hum_interp
                      wvmr_interp = interp_sounding(wvmr_s, pressures_s,y_points)
                      sp_hum_interp = interp_sounding(sp_hum_s, pressures_s,y_points)
                      sat_temp_interp = interp_sounding(sat_temp_s, pressures_s,y_points)
                      theta_e_interp = interp_sounding(theta_e_s, pressures_s,y_points)
                
                      theta_e_sat_interp = interp_sounding(theta_e_sat_s, pressures_s, y_points)
                
                      #virtual_pot_temp_interp = interp_sounding(virtual_pot_temp_s, pressures_s, y_points)
  

                      # Calculate lifting condensation level 
    
   
                      mean_temp_k = MeanFirst500m(temp_interp, geop_height_interp, st_height)
                      dew_mean_temp_k = MeanFirst500m(temp_interp - dewpoints_interp, geop_height_interp, st_height)
                      mean_pres = MeanFirst500m(y_points, geop_height_interp, st_height)

                      lcl_temp = LiftedCondensationLevelTemp(mean_temp_k, dew_mean_temp_k)

                      #wvmr_parcel = WaterVapourMixingRatio(mean_pres, dew_mean_temp_k)
                      #kappa=PoissonConstant(wvmr_parcel)

                      #lcl_pres = LiftedCondensationLevelPres(mean_pres, lcl_temp, mean_temp_k)
    
                      temp_lcl.append(lcl_temp)
                      # Determine planetary boundary layer pressure from parcel vpt method
        
                      #pdb.set_trace()

                      if surface_p_bool == True:

                           pbl_pressure = PBLFromParcelVPT(y_points, temp_interp-273.15, temp_interp - dewpoints_interp-273.15, surface_p, 500.)
               
                           surface_p_l.append(surface_p)

                           surface_p_bool=False
                      else:

                           surface_p_l.append(np.nan)
                           pbl_pressure = np.nan
                           #print pbl_pressure
            
                      pbl_pressure_l.append(pbl_pressure)  
            
          # Save surface pressure and date 
          #surface_p_l.append(surface_p)
           
                      dates_for_plotting_single_single.append(date)
            
                      if len(pbl_pressure_l)>len(temp_lcl):
                           print date
                           print pbl_pressure_l
                           print temp_lcl
                
                      try:
              
                           equil_level_s, parc_prof, lclp_s,lfcp_s, lclt_s, delta_z, CAPE_s, CIN_s = CapeCinPBLInput(
                                                        y_points/100, 
                                                        temp_interp, 
                                                        temp_interp-dewpoints_interp, 
                                                        geop_height_interp, st_height, pbl_pressure/100)  
                
                
                #print CAPE
            
                      except Exception,e:
                           #print e
                           CAPE_s = np.nan
                           CIN_s = np.nan
                           lclp_s = np.nan
                           lclt_s = np.nan
                           lfcp_s = np.nan
                           equil_level_s = np.nan
                           #pass
            
                      CAPE.append(CAPE_s)
                      CIN.append(CIN_s)
                      lclp.append(lclp_s)
                      lclt.append(lclt_s)
                      lfcp.append(lfcp_s)
                      equil_level.append(equil_level_s)
            
                      #Find T=0 pressure
    
                      interp = scipy.interpolate.interp1d(temp_interp ,y_points , bounds_error=False, fill_value=np.nan)
                      t_zero_line.append(float(interp(273.15)))  # Interpolate to T=0
                
 
                    # Find point closest to (surface pressure -50 hPa), if it is not too far away (e.g in a sounding where 
                    # there are no interpolated values near the surface becuase of no data in the sounding near the surface)
                    #  This is forcalculating theta-e saturated minus boundary layer theta-e
    
                      if any(temp_interp[~np.isnan(temp_interp)]) and any(dewpoints_interp[~np.isnan(dewpoints_interp)]) and ('surface_p' in globals()): #If arrays are not totally nan

                         nan_mask = np.ma.masked_array(np.array(y_points), np.isnan(np.array(dewpoints_interp, dtype=float)) | np.isnan(np.array(temp_interp, dtype=float)))

                         i_idx_50hPa = SurfaceMinus(nan_mask, surface_p, 5000.)   

                         if ~np.isnan(i_idx_50hPa):  # If there is a point not too far away from (surface pressure - 50 hPa)
            
                              #Theta-e BL minus theta-es local
    
                              theta_es_min_theta_e_s = theta_e_interp [i_idx_50hPa] - theta_e_sat_interp
                
                              #if np.isnan(sat_temp_interp[i_idx_50hPa]):
                              #    print sat_temp_interp
                              #   st=sat_temp_interp
                              #  vp=vap_press_interp
                              # i_d=i_idx
                              # nm=nan_mask
                              # sf=surf_level
                              # sp=surface_p
                
                              # Makes the next stage simpler if all arrays have the same number of soundings
                
                         else:
                
                              theta_es_min_theta_e_s = np.empty((theta_e_sat_interp.shape))
                              theta_es_min_theta_e_s[:] = np.NAN
                      
            
                         #else:
                         # temp_lcl.append(float('NaN'))
                         # lcl_method_2.append(float('NaN'))
                         # pbl_pressure.append(float('NaN'))
          
                         #else:
                         # temp_lcl.append(float('NaN'))
                         #lcl_method_2.append(float('NaN'))
                         # pbl_pressure.append(float('NaN'))
           
                         try:
                              del surface_p
                         except NameError:
                              pass
                
                         pressures.append(list(y_points))
                         temps.append(list(temp_interp))
                         dewpoints.append(list(dewpoints_interp))
                         winddirs.append(list(winddirs_interp))
                         windspeeds.append(list(windspeeds_interp))
                         geop_height.append(list(geop_height_interp))
                 
                         pot_temp.append(list(pot_temp_interp))
                         sat_vap_pres.append(list(sat_vap_pres_interp))
                         vap_press.append(list(vap_press_interp))
                         rel_hum.append(list(rel_hum_interp))
                         wvmr.append(list(wvmr_interp))
                         sp_hum.append(list(sp_hum_interp))
                         sat_temp.append(list(sat_temp_interp))
                         theta_e.append(list(theta_e_interp))
                
                         theta_e_sat.append(list(theta_e_sat_interp))
                
                         theta_es_min_theta_e.append(list(theta_es_min_theta_e_s))
        
                         dates_for_plotting_single.append(date)
            
                         #virtual_pot_temp.append(list(virtual_pot_temp_interp))
            
                         # dates_for_plotting.append(list([dates_for_plotting_single for i in range(len(y_points))]))
          
                         #if any(dewpoints_interp[~np.isnan(dewpoints_interp)]):
                
                         #pp=pressures_s
                         #tt=temps_s
                         #ppp=pressures
                         #ttt=temps
                         #ti=temp_interp
                         #vpt=virtual_pot_temp_interp
              
                         pressures_s=[]
                         temps_s=[]
                         dewpoints_s=[]
                         winddirs_s=[]
                         windspeeds_s=[]
                         geop_height_s=[]
            
                         temps_s_mask=[]
                         pressures_s_mask =[]
                         dewpoints_s_mask =[]
            
                         temps_scent = []
                         temps_s = []
                         pressures_s_hPa = []
                         dewpoints_s = []
            
                         dewpoint_temp = []
        
                         pot_temp_s = []
    
                         sat_vap_pres_s = []
         
                         vap_press_s = []
            
                         rel_hum_s = []
    
                         wvmr_s =[]
            
                         sp_hum_s = []
            
                         sat_temp_s =[]
        
                         theta_e_s = []
            
                         theta_e_sat_s = []
            
                         theta_es_min_theta_e_s = []
            
                         #virtual_pot_temp_s=[]
            
                         #except ValueError, e:
                         #print e
                         #print line
                         #pass
                         # np.append(pressures_for_plotting, (y_points.tolist(), temp_interp.tolist(), dewpoints_interp.tolist(), winddirs_interp.tolist(), windspeeds_interp.tolist()), axis=1)
                         # levs=0

      
                         #print pressures_s
        

        
                         #print pressures.shape  
             
      #pdb.set_trace()

      u_wind,v_wind = UVWinds(winddirs, windspeeds)
            
      pressures_for_plotting = np.array((pressures, temps, dewpoints, winddirs, windspeeds, pot_temp, sat_vap_pres, vap_press, rel_hum, wvmr, sp_hum, sat_temp, theta_e, theta_e_sat, theta_es_min_theta_e, geop_height, u_wind, v_wind))
      single_value_vars = np.array((temp_lcl, pbl_pressure_l, surface_p_l, t_zero_line, CAPE, CIN, lclp, lclt, lfcp, equil_level))
        
      print pressures_for_plotting.shape
      np.savez('/nfs/a90/eepdw/Data/Observations/Radiosonde_Numpy/Radiosonde_Single_Station_PRESSURES__IND_INTERP_SOUNDING_%s'\
              % (stat), pressures_for_plotting=pressures_for_plotting, single_value_vars=single_value_vars, dates_for_plotting_single=dates_for_plotting_single, dates_for_plotting_single_single=dates_for_plotting_single_single, variable_list=variable_list, variable_list_line=variable_list_line)
     
