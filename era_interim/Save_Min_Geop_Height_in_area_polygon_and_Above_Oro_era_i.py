import cPickle as pickle

from shapely.geometry import Point, Polygon

import numpy as np

from netCDF4 import Dataset

from scipy.interpolate import griddata

p_lev = 925

geopotential, longitude_dom, latitude_dom, time_dom, time_hour  = pickle.load\
                                             (open('/nfs/a90/eepdw/Data/Saved_data/era_i/era_i_emb_time_update_large_geopotential.p', 'rb'))

pressure_levels =  pickle.load(open('/nfs/a90/eepdw/Data/Saved_data/era_i/era_i_emb_pressure_levels.p', 'rb'))

# Monsoon Trough and Ganga Basin combined
polygon = Polygon(((73., 21.), (83., 16.), (87., 22.), (90.,22.), (90.,23.8), (83., 24.2), (76.3, 28.))) 

# Ganga Basin
#polygon = Polygon(((87., 22), (75., 27), (76.3, 30.), (83, 26.2), (90, 25.8), (90., 22)))

lons_data= longitude_dom[0]
lats_data = latitude_dom[0]
lons_data,lats_data = np.meshgrid(lons_data, lats_data)

#  Find points that are within defined polygon

points = np.array([[long,lat] for long, lat in zip(lons_data.flatten(), lats_data.flatten())])
intersects = np.array(map(polygon.intersects, map(Point, points))).reshape(lons_data.shape)

p_lev_idx = np.where(pressure_levels==p_lev)
geopotential_polygon = geopotential[:, p_lev_idx, intersects]


# Do the same for surface geopotential  (still in netcdf format)

#nc = Dataset('/nfs/a90/eepdw/Data/Era_Interim/LandSeaMask/Land_Sea_Mask.nc')
nc = Dataset('/nfs/a90/eepdw/Data/Era_Interim/Orography/era_i_geopotential.nc')

'''
ECMWF give orography as geopotential, which is apparently converted to height by using the WMO gravity constant
9.80665.  Ithought latitude would affect it as well but no mention
http://www.ecmwf.int/en/geopotential-defined-units-m2/s2-both-pressure-levels-and-surface-orography-how-can-height-metres
'''

lons,lats = np.meshgrid(nc.variables['longitude'][:], nc.variables['latitude'][:])
oro_regrid = griddata((lats.flatten(), lons.flatten()), nc.variables['z'][:].flatten(), (lats_data,lons_data), method='linear')

oro_polygon = oro_regrid[intersects]

vals = np.where(geopotential_polygon>(oro_polygon/9.80665), geopotential_polygon, np.nan)

min_geop_full_time = np.min(vals,axis=-1)[:,0]
day_mean_min_geop = [np.mean(min_geop_full_time[np.where(time_dom==day)]) for day in np.unique(time_dom)]

np.savez(
    '/nfs/a90/eepdw/Data/Era_Interim/Era_interim_TimeVar_on_p_levs_mean_by_day_land_domain_'\
    'constrain__and_oro_not_greater_than_data_monsoon_trough_%s' % p_lev, 
    data=day_mean_min_geop, time_coords=np.unique(time_dom), pressures=pressure_levels[p_lev_idx])




