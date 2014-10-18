import numpy as np

from netCDF4 import Dataset
from scipy.interpolate import griddata

from collections import defaultdict

from shapely.geometry import Point, Polygon

from datetime import datetime

lon_high = 101.866 
lon_low = 64.115

lat_high= 33.
lat_low=-6.79

polygon = Polygon(((73., 21.), (83., 16.), (87., 22.), (75., 27.)))

numpy_cube=np.load('/nfs/a90/eepdw/Data/Observations/Satellite/GSMAP_Aug_Sep_2011/GSMAP_EMBRACE.npz')

# Load land sea mask.  TRMM land sea mask is in % of water coverage so 100% is all water

nc = Dataset('/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/TMPA_mask.nc')

#  Regrid lsm to data grid (offset b 0.125 degrees

lsm_lons, lsm_lats = np.meshgrid(nc.variables['lon'][:],nc.variables['lat'][:])
lons_data, lats_data = np.meshgrid(numpy_cube['lons'], numpy_cube['lats'])
#lsm_regrid = griddata((lsm_lats.flatten(), lsm_lons.flatten()), nc.variables['landseamask'][:].flatten(), (lats_data,lons_data), method='linear')
lsm_regrid = griddata((lsm_lats.flatten(), lsm_lons.flatten()), nc.variables['landseamask'][:].flatten(), (lats_data,lons_data), method='linear')

#points = np.array([[lat,lon] for lat, lon in zip(lats_data.flatten(), lons_data.flatten())])
points = np.array([[lon,lat] for lon, lat in zip(lons_data.flatten(), lats_data.flatten())])
intersects = np.array(map(polygon.intersects, map(Point, points))).reshape(lons_data.shape)
pcp_dom_2 = numpy_cube['rain_rates'][intersects,:] #  Need to vary between TRMM, CMORPH etc depending on shape of array
lsm = lsm_regrid[intersects]

bad_values=np.ma.masked_array(pcp_dom_2,pcp_dom_2<0.)

print numpy_cube['rain_rates'].shape
print pcp_dom_2.shape
print lsm.shape

####################################################

# Calculate mean for every time in the date range for entire area

#mean_of_each_time = pcp_dom_2.mean(axis=1)
#mean_and_hour=zip(mean_of_each_time,time_hour)

# OCEAN - Calculate mean for every time in the date range

#lsm_weights=lsm/100

#mean_oc = np.ma.average(pcp_dom_2, axis=1, weights=lsm_weights)
#oc_mean_and_hour=zip(mean_oc,time_hour)

# LAND - Calculate mean for every time in the date range

lsm_weights=1-(lsm/100)
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
np.savez("/nfs/a90/eepdw/Data/Observations/Satellite/GSMAP_Aug_Sep_2011/Diurnal/land_gsmap_diurnal_average_monsoon_trough_polygon" , mean=mean_l, hour=hour_l )
