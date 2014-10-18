import numpy as np
import re
import gzip

from datetime import datetime, timedelta
#from dateutil.relativedelta import relativedelta
from time import mktime

lon_max = 116 
lon_min = 30.5

lat_max= 40
lat_min=-11.25

# Calculate lon/lat grid points from info on website
lons=np.arange(0.05, 360, 0.1)
lons=np.where(lons>180, lons-360, lons)

lats = np.arange(59.95, -59.95-0.1, -0.1)

data_lon_size=3600
data_lat_size=1200

# Get min and max  index positions for latitude and longitude - FOR PCP

la_index_pcp = np.where((lats<=lat_max) & (lats >= lat_min))
lo_index_pcp = np.where((lons<=lon_max) & (lons >= lon_min))

la_i_max_pcp = np.max(la_index_pcp)
la_i_min_pcp = np.min(la_index_pcp)
lo_i_max_pcp = np.max(lo_index_pcp)
lo_i_min_pcp = np.min(lo_index_pcp)


print la_i_min_pcp,la_i_max_pcp, lo_i_min_pcp,lo_i_max_pcp

date_max=datetime(2011, 9, 8, 0, 0)
date_min=datetime(2011, 8, 18, 0, 0)

date_max_unix=mktime(date_max.timetuple())
date_min_unix=mktime(date_min.timetuple())
number_of_hours=(date_max_unix-date_min_unix)/3600

no_of_lon_coords=lo_i_max_pcp-lo_i_min_pcp+1
no_of_lat_coords=la_i_max_pcp-la_i_min_pcp+1

dt = np.dtype('<f4') # Little-endian 4byte (32bit) float

rain_rates=np.empty((no_of_lat_coords,no_of_lon_coords,number_of_hours+1), dtype=dt) # Data
#rain_rates_full=np.empty((3600,1200,number_of_hours+1), dtype=dt) # Data
#rain_rates_full=np.empty((1200,3600,number_of_hours+1), dtype=dt) # Data
time_list=np.empty((4, number_of_hours+1)) 

data_dir = '/nfs/a90/eepdw/Data/Observations/Satellite/GSMAP_Aug_Sep_2011/realtime/archive'

date_t=date_min

time_count=0
while date_t<=date_max:
    print 'Processing %s/%s/%s/%s/gsmap_nrt.%s%s%s.%s.dat.gz' % (data_dir, 
                  date_t.year, date_t.strftime('%m'), date_t.strftime('%d'),
                  date_t.year, date_t.strftime('%m'), date_t.strftime('%d'), date_t.strftime('%H%M'))
    f = gzip.open('%s/%s/%s/%s/gsmap_nrt.%s%s%s.%s.dat.gz' % (data_dir, 
                  date_t.year, date_t.strftime('%m'), date_t.strftime('%d'),
                  date_t.year, date_t.strftime('%m'), date_t.strftime('%d'), date_t.strftime('%H%M')), 'rb')
    #rain_rates[:,:,time_count] = np.frombuffer(f.read(), dtype=dt).reshape(data_lon_size,data_lat_size)[lo_i_min_pcp:lo_i_max_pcp+1, la_i_min_pcp:la_i_max_pcp+1]
    rain_rates[:,:,time_count] = np.frombuffer(f.read(), dtype=dt).reshape(data_lat_size,data_lon_size)[la_i_min_pcp:la_i_max_pcp+1, lo_i_min_pcp:lo_i_max_pcp+1]
    #rain_rates_full[:,:,time_count] = np.frombuffer(f.read(), dtype=dt).reshape(data_lon_size,data_lat_size)
    #rain_rates_full[:,:,time_count] = np.frombuffer(f.read(), dtype=dt).reshape(data_lat_size,data_lon_size)
    f.close()

    time_list[:, time_count] = (date_t.year, date_t.month, date_t.day, date_t.hour)  
    
    date_t=date_t+timedelta(hours=1)
    time_count+=1

np.savez('/nfs/a90/eepdw/Data/Observations/Satellite/GSMAP_Aug_Sep_2011/GSMAP_EMBRACE', rain_rates=rain_rates, time_list=time_list, lats=lats[la_index_pcp], lons=lons[lo_index_pcp])
#np.savez('/nfs/a90/eepdw/Data/Observations/Satellite/GSMAP_Aug_Sep_2011/GSMAP_EMBRACE_full', rain_rates=rain_rates_full, time_list=time_list, lats=lats[la_index_pcp], lons=lons[lo_index_pcp])
