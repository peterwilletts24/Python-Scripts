# Plot average TRMM precipitation for embrace period on map
#

def TRMM_Average():
 import matplotlib.pyplot as plt
 import numpy as np
 from netCDF4 import Dataset
 import os

# constants
# lon from 60 105
# lat -10 30 

 lon_amounts = 180
 lat_amounts = 160
 dys = 21
 yy = 2011
 #hrs = ['00','03','06','09','12','15','18','21']
 hrs = [0,3,6,9,12,15,18,21]
 dr = '/nfs/a80/earceb/model_data/observations/satellite/TRMM/2011/'
 
#create master array
 diurn = np.zeros((8),dtype = float)
 mastr = np.zeros((lat_amounts, lon_amounts, dys),dtype =int)
 
 mnth= 8
 for dte in['18','19','20','21','22','23','24','25','26','27','28','29','30','31']:
#will first test if file exists
  for h in hrs:
   fle = str(dr)+'3B42.'+str(yy)+'0'+str(mnth)+str(dte)+'.'+str(h)+'.7.nc'
   #print fle
   if os.path.isfile(fle):
    print fle 
    # create holder array for day date
    prec_dat = np.zeros((lat_amounts,lon_amounts,8),dtype = float)
    nc = Dataset(str(dr)+'3B42.'+str(yy)+'0'+str(mnth)+str(dte)+'.'+str(h)+'.7.nc')
    trmm = nc.variables['pcp']
    print trmm
    mat = trmm[1, 160:320, 960:1140]
    nc.close()

    # transfer data into prec_dat
    #prec_dat[:,:,h/3] = mat[:]
    #for lat in range(0,lat_amounts):
     #for lon in range(0,lon_amounts):
      #if prec_dat[lat,lon,h/3] < 0:
       #continue
      #else:
       # diurn[h/3] = diurn [h/3] + prec_dat[lat,lon,h/3]
     #print diurn[h/3]
   #else:
   # print 'Error'+str(dte)
