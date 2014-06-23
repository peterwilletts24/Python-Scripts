########################################
# Check station metadata text file to find stations within lat/lon
#
# Created by: Peter Willetts
# Created on: 14/5/2014
#
########################################
#
# http://www1.ncdc.noaa.gov/pub/data/igra/readme.txt
###################################################

import glob
import re
import numpy as np

lon_max = 116 
lon_min = 30.5

lat_max= 40
lat_min=-11.25

#lon_max = 102 
#lon_min = 64

#lat_max= 30
#lat_min=-10.5

station_list='/nfs/a90/eepdw/Data/Observations/Radiosonde_downloaded_from_NOAA_GUAN/igra-stations.txt'

station_metadata=[]
f = open(station_list,'r')

station_list_inside_domain=[]
for line in f:
     #line2 = [s.strip() for s in line.split('  ') if s.strip()]
     line=line.strip()
     line=re.sub(r'([A-Z])\s([A-Z])', r'\1_\2',line)
     line=re.sub(r'([A-Z])\s\s([A-Z])', r'\1_\2',line)
     line=line.split()
     try:         
      #float(line[-3])
      #print line
      #print (line[4])
      if (float(line[3])<lat_max and float(line[3])>lat_min and float(line[4])<lon_max and float(line[4])>lon_min and float(line[-1])>=2011):
             station_list_inside_domain.append(line)
     except ValueError:
         #try:
         #    float(line[-4])
          #   if (float(line[-6])<lat_max and float(line[-6])>lat_min and float(line[-5])<lon_max and float(line[-5])>lon_min and float(line[-1])>=2011):
         #        station_list_inside_domain.append(line)  
         #except ValueError:
          #   try:
            #     if (float(line[-7])<lat_max and float(line[-7])>lat_min and float(line[-6])<lon_max and float(line[-6])>lon_min and float(line[-1])>=2011):
            #         station_list_inside_domain.append(line)
            # except ValueError:
            #      station_list_inside_domain.append('Failed lat/lon check - %s' % line)
         print line[-5]
         print len(line)
         
print station_list_inside_domain

dataFile = open('/nfs/see-fs-01_users/eepdw/python_scripts/radiosonde_data/stations_inside_lat_lon_box_large.txt', 'w')
#dataFile = open('/nfs/see-fs-01_users/eepdw/python_scripts/radiosonde_data/stations_inside_lat_lon_box.txt', 'w')
for eachitem in station_list_inside_domain:
 #   print eachitem
    dataFile.write(str(eachitem)+'\n')

     #print columns[4]
     #print columns

dataFile.close()
