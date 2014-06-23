import cPickle as pickle
import numpy as np
from collections import defaultdict
from netCDF4 import Dataset
from scipy.interpolate import griddata

from shapely.geometry import Point, Polygon

polygon = Polygon(((77, 10), (76., 16.), (80., 16.), (80., 10.)))
pcp_dom, longitude_dom, latitude_dom, time_dom, time_hour = pickle.load(open('/nfs/a90/eepdw/Data/Saved_data/TRMM/trmm_emb_time_update_large.p', 'rb'))

# Load land sea mask.  TRMM land sea mask is in % of water coverage so 100% is all water

nc = Dataset('/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/TMPA_mask.nc')

#  Regrid lsm to data grid (offset b 0.125 degrees

lsm_lons, lsm_lats = np.meshgrid(nc.variables['lon'][:],nc.variables['lat'][:])
lons_data, lats_data = np.meshgrid(longitude_dom[0], latitude_dom[0])
lsm_regrid = griddata((lsm_lats.flatten(), lsm_lons.flatten()), nc.variables['landseamask'][:].flatten(), (lats_data,lons_data), method='linear')

#  Find points that are within defined polygon

points = np.array([[long,lat] for long, lat in zip(lons_data.flatten(), lats_data.flatten())])
intersects = np.array(map(polygon.intersects, map(Point, points))).reshape(lons_data.shape)
pcp_dom_2 = pcp_dom[:,intersects]
lsm = lsm_regrid[intersects]

print pcp_dom.shape
print pcp_dom_2.shape
print lsm.shape
####################################################

# Calculate mean for every time in the date range for entire area

mean_of_each_time = pcp_dom_2.mean(axis=1)
mean_and_hour=zip(mean_of_each_time,time_hour)

# OCEAN - Calculate mean for every time in the date range

lsm_weights=lsm/100

mean_oc = np.ma.average(pcp_dom_2, axis=1, weights=lsm_weights)
oc_mean_and_hour=zip(mean_oc,time_hour)

# LAND - Calculate mean for every time in the date range

lsm_weights=1-(lsm/100)
mean_la = np.ma.average(pcp_dom_2, weights=lsm_weights, axis=1)
la_mean_and_hour=zip(mean_la,time_hour)

#####################################################


# Sort into time of day #################

# Total

i = defaultdict(list)

for v,k in mean_and_hour:
    i[k.strip()].append(v)

mean=[]
hour=[]

# Average for each time of day

for q,a in enumerate(i.items()):
    print a[1]
    if a[1]:
        mean.append(np.mean(a[1]))
        hour.append(a[0])

print mean
print hour

# Land

i = defaultdict(list)

for v,k in la_mean_and_hour:
    i[k.strip()].append(v)

mean_l=[]
hour_l=[]

# Average for each time of day

for q,a in enumerate(i.items()):
    print a[1]
    if a[1]:
        mean_l.append(np.mean(a[1]))
        hour_l.append(a[0])

print mean_l
print hour_l

# Ocean

i = defaultdict(list)

for v,k in oc_mean_and_hour:
    i[k.strip()].append(v)

mean_o=[]
hour_o=[]

      # Average for each time of day

for q,a in enumerate(i.items()):
    #print a[1]
    if a[1]:
        mean_o.append(np.mean(a[1]))
        hour_o.append(a[0])

print mean_o
print hour_o

# Save

np.savez("/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/Diurnal/total_trmm_diurnal_average_rain_shadow_polygon", mean=mean, hour=hour )
#np.savez("/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/Diurnal/sea_trmm_diurnal_average_rain_shadow_polygon" , mean=mean_o, hour=hour_o )
np.savez("/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/Diurnal/land_trmm_diurnal_average_rain_shadow_polygon" , mean=mean_l, hour=hour_l )
