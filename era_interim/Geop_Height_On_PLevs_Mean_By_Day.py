import cPickle as pickle

from shapely.geometry import Point, Polygon

import numpy as np

from netCDF4 import Dataset

from scipy.interpolate import griddata

import pdb


geopotential, longitude_dom, latitude_dom, time_dom, time_hour  = pickle.load\
                                             (open('/nfs/a90/eepdw/Data/Saved_data/era_i/era_i_emb_time_update_large_geopotential.p', 'rb'))

pressure_levels =  pickle.load(open('/nfs/a90/eepdw/Data/Saved_data/era_i/era_i_emb_pressure_levels.p', 'rb'))

lons_data= longitude_dom[0]
lats_data = latitude_dom[0]
lons_data,lats_data = np.meshgrid(lons_data, lats_data)

#pdb.set_trace()

day_mean_min_geop = [np.mean(geopotential[np.where(time_dom==day)], axis=0) for day in np.unique(time_dom)]

np.savez(
    '/nfs/a90/eepdw/Data/Era_Interim/Era_Interim_Daily_Geopotential_Height_EMBRACE_Period'\
    , data=day_mean_min_geop, time_coords=np.unique(time_dom), pressures=pressure_levels, longitudes = lons_data, latitudes = lats_data)


