import numpy as np

from netCDF4 import Dataset
from scipy.interpolate import griddata

from collections import defaultdict

from shapely.geometry import Point, Polygon

import datetime

import iris
import cPickle as pickle


import pdb

lon_high = 101.866 
lon_low = 64.115

lat_high= 33.
lat_low=-6.79

polygon = Polygon(((73., 21.), (83., 16.), (87., 22.), (75., 27.)))

numpy_cube=np.load('/nfs/a90/eepdw/Data/Observations/Satellite/GSMAP_Aug_Sep_2011/GSMAP_EMBRACE.npz')

trmm_pcp, longitude_trmm, latitude_trmm, time_trmm, time_hour_trmm = pickle.load(open('/nfs/a90/eepdw/Data/Saved_data/TRMM/trmm_emb_time_update_large.p', 'rb'))
 
hours_in_trmm = np.unique(np.array(time_hour_trmm, dtype=int))
lats_in_trmm = np.unique(latitude_trmm)
lons_in_trmm = np.unique(longitude_trmm)


# Load land sea mask.  TRMM land sea mask is in % of water coverage so 100% is all water

#nc = Dataset('/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/TMPA_mask.nc')

lsm_cube =iris.load_cube('/nfs/a90/eepdw/Data/EMBRACE/dkbh/dkbhu/30.pp', 'land_binary_mask')

#  Regrid lsm to data grid (offset b 0.125 degrees

#lsm_lons, lsm_lats = np.meshgrid(nc.variables['lon'][:],nc.variables['lat'][:])

lat = lsm_cube.coord('grid_latitude').points
lon = lsm_cube.coord('grid_longitude').points
                    
cs = lsm_cube.coord_system('CoordSystem')
                    
lons, lats = np.meshgrid(lon, lat) 
lsm_lons, lsm_lats = iris.analysis.cartography.unrotate_pole\
                                 (lons,lats, cs.grid_north_pole_longitude, cs.grid_north_pole_latitude)


lons_data, lats_data = np.meshgrid(numpy_cube['lons'], numpy_cube['lats'])
#lsm_regrid = griddata((lsm_lats.flatten(), lsm_lons.flatten()), nc.variables['landseamask'][:].flatten(), (lats_data,lons_data), method='linear')
print lons_data.shape
lsm_regrid = griddata((lsm_lats.flatten(), lsm_lons.flatten()), lsm_cube.data.flatten(), (lats_data,lons_data), method='linear')

#lsm_regrid = griddata((lsm_lats.flatten(), lsm_lons.flatten()), nc.variables['landseamask'][:].flatten(), (lats_data,lons_data), method='linear')

#pdb.set_trace()

datetimes_list = [datetime.datetime(int(dt[0]), int(dt[1]), int(dt[2]), int(dt[3])) for dt in np.rollaxis(numpy_cube['time_list'], 1)]
timestamp_list  = [(dt - datetime.datetime(1970, 1, 1)).total_seconds() for dt in datetimes_list]

trmm_datetime_list = [datetime.datetime.strptime(tt, '%Y-%m-%d %H' ) for tt in time_trmm]
trmm_timestamp_list  = [(dt - datetime.datetime(1970, 1, 1)).total_seconds() for dt in trmm_datetime_list]

pdb.set_trace()

#points = np.array([[lat,lon] for lat, lon in zip(lats_data.flatten(), lons_data.flatten())])

bad_values=np.ma.masked_array(numpy_cube['rain_rates'], numpy_cube['rain_rates']<0.)

data_regridded_to_trmm = griddata((np.tile(lats_data.flatten(), len(timestamp_list)), np.tile(lons_data.flatten(), len(timestamp_list)), 
                            np.tile(timestamp_list, lons_data.flatten().shape)), 
                            bad_values.flatten(), (longitude_trmm, latitude_trmm, trmm_timestamp_list), method='linear')

points = np.array([[lon,lat] for lon, lat in zip(longitude_trmm.flatten(), latitude_trmm.flatten())])
intersects = np.array(map(polygon.intersects, map(Point, points))).reshape(lons_data.shape)
pcp_dom_2 = numpy_cube['rain_rates'][intersects,:] #  Need to vary between TRMM, CMORPH etc depending on shape of array
lsm = lsm_regrid[intersects]

print numpy_cube['rain_rates'].shape
print pcp_dom_2.shape
print lsm.shape

#pdb.set_trace()


####################################################

# Calculate mean for every time in the date range for entire area

#mean_of_each_time = pcp_dom_2.mean(axis=1)
#mean_and_hour=zip(mean_of_each_time,time_hour)

# OCEAN - Calculate mean for every time in the date range

#lsm_weights=lsm/100

#mean_oc = np.ma.average(pcp_dom_2, axis=1, weights=lsm_weights)
#oc_mean_and_hour=zip(mean_oc,time_hour)

# LAND - Calculate mean for every time in the date range
#lsm_weights=1-(lsm/100) # Using TRMM mask which is in range 0 to 100
#pdb.set_trace()

lsm_weights=lsm  # Using UM LSM which is in range 0 to 1

mean_la = np.ma.average(bad_values, weights=lsm_weights, axis=0)
#pdb.set_trace()
time_hour =  numpy_cube['time_list'][3]
#pdb.set_trace()
la_mean_and_hour=zip(mean_la.data.flatten(),time_hour)

#####################################################


# Sort into time of day #################

# Total

#i = defaultdict(list)

#for v,k in mean_and_hour:
#    i[k.strip()].append(v)

#mean=[]
#hour=[]

# Average for each time of day

#for q,a in enumerate(i.items()):
    #print a[1]
 #   if a[1]:
#        mean.append(np.mean(a[1]))
#        hour.append(a[0])

#print mean
#print hour

# Land

i = defaultdict(list)

for v,k in la_mean_and_hour:
    i[k].append(v)

mean_l=[]
hour_l=[]

# Average for each time of day

for q,a in enumerate(i.items()):
    #print a[1]
    if a[1]:
        mean_l.append(np.mean(a[1]))
        hour_l.append(a[0])

print mean_l
print hour_l

# Ocean

#i = defaultdict(list)

#for v,k in oc_mean_and_hour:
 #   i[k.strip()].append(v)

#mean_o=[]
#hour_o=[]

      # Average for each time of day

#for q,a in enumerate(i.items()):
    #print a[1]
#    if a[1]:
#        mean_o.append(np.mean(a[1]))
#        hour_o.append(a[0])

#print mean_o
#print hour_o

# Save

#np.savez("/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/Diurnal/total_trmm_diurnal_average_lat_%s_%s_lon_%s_%s_southern_western_indian_ocean" % (lat_low,lat_high, lon_low, lon_high), mean=mean, hour=hour )
#np.savez("/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/Diurnal/sea_trmm_diurnal_average_lat_%s_%s_lon_%s_%s_southern_western_indian_ocean" % (lat_low,lat_high, lon_low, lon_high), mean=mean_o, hour=hour_o )
np.savez("/nfs/a90/eepdw/Data/Observations/Satellite/GSMAP_Aug_Sep_2011/Diurnal/land_gsmap_diurnal_average_monsoon_trough_polygon_3hourly_point25deg" , mean=mean_l, hour=hour_l )






