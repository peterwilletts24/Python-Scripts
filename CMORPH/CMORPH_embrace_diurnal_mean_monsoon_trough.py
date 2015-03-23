import cPickle as pickle
import numpy as np
from collections import defaultdict
from netCDF4 import Dataset
from scipy.interpolate import griddata

from shapely.geometry import Point, Polygon

from datetime import datetime

import pdb

import iris

# CMOPRH is 0 to 360 longitude . . .

polygon = Polygon(((73., 21.), (83., 16.), (87., 22.), (75., 27.)))

pcp_dom, longitude_dom, latitude_dom, time_dom = pickle.load(open('/nfs/a90/eepdw/Data/Saved_data/CMORPH/cmorph_emb_time_update_large.p', 'rb'))

# Load land sea mask.  TRMM land sea mask is in % of water coverage so 100% is all water
#
#
#
#nc = Dataset('/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/TMPA_mask.nc')

#lsm_lons, lsm_lats = np.meshgrid(nc.variables['lon'][:],nc.variables['lat'][:])

lsm_cube =iris.load_cube('/nfs/a90/eepdw/Data/EMBRACE/dkbh/dkbhu/30.pp', 'land_binary_mask')

def LSM_Regrid_From_IRIS_Cube(lsm_cube, lats_data, lons_data):   
    '''
    Regrid lsm iris cube to numpy array data grid
    '''

    lat = lsm_cube.coord('grid_latitude').points
    lon = lsm_cube.coord('grid_longitude').points
                    
    cs = lsm_cube.coord_system('CoordSystem')
                    
    lons, lats = np.meshgrid(lon, lat) 
    lsm_lons, lsm_lats = iris.analysis.cartography.unrotate_pole\
                                 (lons,lats, cs.grid_north_pole_longitude, cs.grid_north_pole_latitude)

    lsm_lons = np.where(lsm_lons>180, lsm_lons-360, lsm_lons)

    lsm_regrid = griddata((lsm_lats.flatten(), lsm_lons.flatten()), lsm_cube.data.flatten(), (lats_data,lons_data), method='linear')


    return lsm_regrid

lons_data, lats_data = np.meshgrid(longitude_dom[0], latitude_dom[0])

lons_data = np.where(lons_data>180, lons_data-360, lons_data)

    #pdb.set_trace()

print lons_data.shape
    
#pdb.set_trace()

lsm_regrid = LSM_Regrid_From_IRIS_Cube(lsm_cube, lats_data, lons_data)

print lsm_regrid.shape
#pdb.set_trace()
#  Find points that are within defined polygon

#points = np.array([[lat,lon] for lat, lon in zip(lats_data.flatten(), lons_data.flatten())])
points = np.array([[long,lat] for long, lat in zip(lons_data.flatten(), lats_data.flatten())])
intersects = np.array(map(polygon.intersects, map(Point, points))).reshape(lons_data.shape)
pcp_dom_2 = pcp_dom[:,:,intersects] #  Need to vary between TRMM, CMORPH etc depending on shape of array
lsm = lsm_regrid[intersects]

bad_values=np.ma.masked_array(pcp_dom_2,pcp_dom_2<0.)

#pdb.set_trace()

print pcp_dom.shape
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

#lsm_weights=1-(lsm/100) # Using TRMM mask which is in range 0 to 100
lsm_weights=lsm  # Using UM LSM which is in range 0 to 1
mean_la = np.ma.average(bad_values, weights=lsm_weights, axis=-1)
#pdb.set_trace()
time_hour = [datetime.utcfromtimestamp(v).hour for v in time_dom.flatten()]
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
np.savez("/nfs/a90/eepdw/Data/Observations/Satellite/CMORPH/Diurnal/land_cmorph_diurnal_average_monsoon_trough_polygon" , mean=mean_l, hour=hour_l )
