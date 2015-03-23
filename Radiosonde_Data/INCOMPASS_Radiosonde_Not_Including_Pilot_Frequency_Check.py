  #Raw data read for sounding date/time - Climatology
    # /nfs/a90/eepdw/Data/Observations/Radiosonde_downloaded_from_NOAA_GUAN/igra-stations.txt
    # /nfs/a90/eepdw/Data/Observations/Radiosonde_downloaded_from_NOAA_GUAN/readme.txt
    
    # Station numbers for Bombay - 43003; New Delhi - 42182, Lucknow - 42369, Cochin, Bangalore, Madras, Minicoy, Nagpur Sonegaon)
import datetime
import scipy.interpolate
import re
import numpy as np

import pdb
    
import os 
    
   
#station_list_cs=[43371, 43285, 43295, 43296, 43279, 43128, 43003, 42647, 42339, 42369, 42182, 42809, 42492, 42361, 42867, 43333, 42971, 42379, 42492]
station_list_cs=[42379]

def UnixEpochWeeklyGrid(date_min, date_max):

    delta = relativedelta(weeks=+1)

    d = date_min

    grid=[]
    while ((d <= date_max) and (d.timetuple() not in grid)): 
        grid.append(timegm(d.timetuple()))
        d += delta
    
    return grid

import re

station_list_search='/nfs/a90/eepdw/Data/Observations/Radiosonde_downloaded_from_NOAA_GUAN/igra-stations.txt'
station_metadata=[]
f = open(station_list_search,'r')
for line in f:
     line = line.strip()
     line=re.sub(r'([A-Z])\s([A-Z])', r'\1_\2',line)
     line=re.sub(r'([A-Z])\s\s([A-Z])', r'\1_\2',line)
     station_metadata.append(line.split())
f.close()

def station_info_search(stat):
    for line in station_metadata:
         if "%s" % stat in line: 
             st = line[2].lower().title().replace('_',' ')
             lo = float(line[3])
             la = float(line[4])
             st_height = float(line[5])
    return st,la,lo, st_height

test_count=0
    
    #match_header = re.compile(r'(#.....20..0[56789|#.....19..0[56789])')
match_header = re.compile(r'(#.....20|#.....19)')
    #match_header_20th = re.compile(r'#.....19..0[56789]')
match_bad = re.compile(r'9999')
        
lev_count=False
    
 
for stat in station_list_cs:
        
      st,la,lo, st_height = station_info_search(stat)
        
      lev_count=False
        
 
      dates_for_plotting_single=[]
      dates_for_plotting_single_single=[]
        
      levs=-1
      numlevs=0
    
      #print stat
    
      i = '/nfs/a90/eepdw/Data/Observations/Radiosonde_downloaded_from_NOAA_GUAN/%s.dat' % stat
      print i
      f = open(i,'r')
      #print f
      for line in f:
       #line = line.strip().replace('-99999',' %s ' %float('NaN'))
       line = line.strip()
       #print line
       #columns=line.split()
                
       if levs >= numlevs and pressures_s and list(temps_s) and list(dewpoints_s) and winddirs_s and windspeeds_s:
          lev_count=False
           
          if ~np.all(np.isnan(pressures_s)) and ~np.all(np.isnan(temps_s)) and ~np.all(np.isnan(dewpoints_s)) and ~np.all(np.isnan(winddirs_s)) and ~np.all(np.isnan(windspeeds_s)):
              dates_for_plotting_single_single.append(date)
                

          pressures_s=[]
          temps_s=[]
          dewpoints_s=[]
          winddirs_s=[]
          windspeeds_s=[]
          geop_height_s=[]
                
       else:
          pdb.set_trace()
        
   
       if lev_count is True: 

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
      
        #print pressures_s
        
    # Extract independent headers 
    
  
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
        
 
             
      #np.savez('/nfs/a90/eepdw/Data/Observations/Radiosonde_Numpy/Radiosonde_Single_Station_DateTimes_Only_%s_%s_to_%s'\
        #      % (stat, date_min.strftime('%Y%m%d'), date_max.strftime('%Y%m%d') ), dates_for_plotting_single_single=dates_for_plotting_single_single)
     
