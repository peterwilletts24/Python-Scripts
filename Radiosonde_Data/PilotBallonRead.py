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

imp.load_source('GenMeteoFuncs', '/nfs/see-fs-01_users/eepdw/python_scripts/modules/GeneralMeteoFunctions.py')
from GenMeteoFuncs import *
imp.load_source('SoundingRoutines', '/nfs/see-fs-01_users/eepdw/python_scripts/Tephigram/Sounding_Routines.py')
from SoundingRoutines import *
imp.load_source('GeogFuncs', '/nfs/see-fs-01_users/eepdw/python_scripts/modules/GeogFunctions.py')
from GeogFuncs import *

station_list_cs=[43003, 43014, 42867, 43371, 43353, 43285,  43192, 43150, 42339, 40990, 40948]
station_list_cs=[42379, 42492, 42475]

y_points=np.linspace(5000, 108000, 200) # Points for pressure interpolation

#match_header = re.compile(r'(#.....20..0[56789|#.....19..0[56789])')
match_header = re.compile(r'(#.....20|#.....19)')
#match_header_20th = re.compile(r'#.....19..0[56789]')
match_bad = re.compile(r'9999')

variable_list={'pressures': 0, 'winddirs':1, 'windspeeds':2, 'geopheight':3, 'u_wind':4, 'v_wind':5}

    
test_count=0
        
for stat in station_list_cs:
        
      i = '/nfs/a90/eepdw/Data/Observations/Radiosonde_downloaded_from_NOAA_GUAN/%s.dat' % stat

      st,la,lo, st_height = StationInfoSearch(stat)
        
      lev_count=False
        
    
      pressures=[]
      winddirs=[]
      windspeeds=[]
      geop_height=[]
        
  
      dates_for_plotting_single=[]
   
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
                winddirs_s=[]
                windspeeds_s=[]
                geop_height_s=[]

            if (lev_count is True) & (match_header.search(line) is None): 
        
             
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
                    winddirs_s.append(float(line[26:31].replace('-9999','%s' %float('NaN')).replace('-8888','%s' %float('NaN')))) 
                    windspeeds_s.append(float(line[31:36].replace('-9999','%s' %float('NaN')).replace('-8888','%s' %float('NaN'))))
                    geop_height_s.append(float(line[9:14].replace('-9999','%s' %float('NaN')).replace('-8888','%s' %float('NaN'))))
         
       
                    levs+=1
                except Exception,e:
                    print e
                    print line


            if (match_header.search(line) is None) and levs >= numlevs and pressures_s and winddirs_s and windspeeds_s:
                    
                 lev_count=False
        
                 if ~np.all(np.isnan(pressures_s))\
                    and ~np.all(np.isnan(winddirs_s)) and ~np.all(np.isnan(windspeeds_s)):
        
                      # Calculate sounding variables

                      geop_height_s = np.ma.masked_array(np.array(geop_height_s, dtype=float), np.isnan(np.array(geop_height_s, dtype=float)))
                      
                      winddirs_interp= interp_sounding(np.array(winddirs_s, dtype=float), pressures_s, y_points)
                      windspeeds_interp= interp_sounding(np.array(windspeeds_s, dtype=float)/10, pressures_s, y_points)
                  
                      geop_height_interp=interp_sounding(geop_height_s.data[~geop_height_s.mask], np.array(pressures_s)[~geop_height_s.mask], y_points)
                      
               
                      try:
                           del surface_p
                      except NameError:
                           pass
                
                      pressures.append(list(y_points))
                   
                      winddirs.append(list(winddirs_interp))
                      windspeeds.append(list(windspeeds_interp))
                      geop_height.append(list(geop_height_interp))
                     
                      dates_for_plotting_single.append(date)
            
                    
              
                      pressures_s=[]
                      
                      winddirs_s=[]
                      windspeeds_s=[]
                      geop_height_s=[]
            
        
                      #pdb.set_trace()

      u_wind,v_wind = UVWinds(winddirs, windspeeds)
            
      pressures_for_plotting = np.array((pressures, winddirs, windspeeds, geop_height, u_wind, v_wind))
     
      print pressures_for_plotting.shape
      np.savez('/nfs/a90/eepdw/Data/Observations/Radiosonde_Numpy/Radiosonde_Single_Station_WINDDATA__IND_INTERP_SOUNDING_%s'\
              % (stat), pressures_for_plotting=pressures_for_plotting,dates_for_plotting_single=dates_for_plotting_single, variable_list=variable_list)
     
