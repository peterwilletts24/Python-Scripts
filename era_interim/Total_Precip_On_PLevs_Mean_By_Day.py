import numpy as np

from netCDF4 import Dataset

from netCDF4 import num2date

import datetime

nc=Dataset('/nfs/a90/eepdw/Data/Era_Interim/Precipitation/Era_Interim_Precip_data_EMBRACE_Period.nc')

times = nc.variables['time']

times = num2date(times[:],units=times.units,calendar=times.calendar)

date_strings= np.array(['%s-%s-%s' % (t.year, t.month, t.day) for t in times])

unique_dates = np.unique(['%s-%s-%s' % (t.year, t.month, t.day) for t in times])

lat_lon_mean=[]

lons_data,lats_data = np.meshgrid( np.array(nc.variables['longitude'][:]), np.array(nc.variables['latitude'][:]))

for ud in unique_dates:

    lat_lon_mean.append(np.mean(nc.variables['tp'][date_strings==ud], axis=0))

np.savez(    
               '/nfs/a90/eepdw/Data/Era_Interim/Era_Interim_Daily_Total_Precip_EMBRACE_Period'\
              , data=np.array(lat_lon_mean), time_coords=unique_dates, longitudes = lons_data, latitudes = lats_data)


