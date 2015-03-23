import cPickle as pickle
import numpy as np
from collections import defaultdict
from netCDF4 import Dataset
from scipy.interpolate import griddata

from shapely.geometry import Point, Polygon

import datetime

import pdb

import iris


# Min and max lats lons from smallest model domain (dkbhu) - see spreadsheet

# Unrotated - can only be used (without unrotating first) for djznw  - global model

latmin=-6.79
latmax=33.04
lonmin=64.12
lonmax=101.87

dtmindt = datetime.datetime(2011,8,28,0,0,0)
dtmaxdt = datetime.datetime(2011,9,2,23,0,0)

pcp_dom, longitude_dom, latitude_dom, time_dom = pickle.load(open('/nfs/a90/eepdw/Data/Saved_data/CMORPH/cmorph_emb_time_update_large.p', 'rb'))

# Load land sea mask.  TRMM land sea mask is in % of water coverage so 100% is all water

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

#  Regrid lsm to data grid (offset b 0.125 degrees

#lsm_lons, lsm_lats = np.meshgrid(nc.variables['lon'][:],nc.variables['lat'][:])
#lons_data, lats_data = np.meshgrid(longitude_dom[0], latitude_dom[0])
#lsm_regrid = griddata((lsm_lats.flatten(), lsm_lons.flatten()), nc.variables['landseamask'][:].flatten(), (lats_data,lons_data), method='linear')

#  Find points that are within defined area

lat_max_idx = np.max(np.where((latmin<=latitude_dom[0]) & (latitude_dom[0]<=latmax)))
lat_min_idx = np.min(np.where((latmin<=latitude_dom[0]) & (latitude_dom[0]<=latmax)))

lon_max_idx = np.max(np.where((lonmin<=longitude_dom[0]) & (longitude_dom[0]<=lonmax)))
lon_min_idx = np.min(np.where((lonmin<=longitude_dom[0]) & (longitude_dom[0]<=lonmax)))

#pdb.set_trace()

datetime_time_dom_list = np.array([datetime.datetime.utcfromtimestamp(v) for v in time_dom.flatten()])
#pdb.set_trace()

time_max_idx = np.max(np.where((datetime_time_dom_list>=dtmindt) & (datetime_time_dom_list<=dtmaxdt)))
time_min_idx = np.min(np.where((datetime_time_dom_list>=dtmindt) & (datetime_time_dom_list<=dtmaxdt)))

pcp_dom = pcp_dom.reshape((datetime_time_dom_list.shape[0], pcp_dom.shape[2], pcp_dom.shape[3]))

pcp_dom_2 = pcp_dom[time_min_idx:time_max_idx, lat_min_idx:lat_max_idx, lon_min_idx:lon_max_idx]
lsm = lsm_regrid[lat_min_idx:lat_max_idx, lon_min_idx:lon_max_idx]

bad_values=np.ma.masked_array(pcp_dom_2,pcp_dom_2<0.)

print pcp_dom.shape
print pcp_dom_2.shape
print lsm.shape
lsm=np.resize(lsm, pcp_dom_2.shape)[0]

####################################################

bins=np.linspace(0.,200., 2000.)

for ls in ['land', 'sea']:
    if ls=='land':
        pdf=[]
        bmin=0.
        for b in bins:
            pdf.append(((np.resize(np.abs(100-lsm), bad_values.shape)[(bad_values <= b) & (bad_values>=bmin)])/100).sum())
            if bmin==0.:
                bmin=0.0000000001
            else:
                bmin=b

    elif ls=='sea':
        pdf=[]
        bmin=0.
        for b in bins:
            pdf.append(((np.resize(lsm, bad_values.shape)[(bad_values <= b) & (bad_values >= bmin)])/100).sum())
            if bmin==0.:
                bmin=0.0000000001
            else:
                bmin=b
#####################################################


# Save


    np.savez("/nfs/a90/eepdw/Data/Observations/Satellite/CMORPH/Rainfall_Intensity/CMORPH_EMBRACE_Rainfall_Intensity_PDF_dkbhu_latlon_%s_%s_%s" 
               % (ls,  dtmindt.strftime('%Y%b%d%H%M'), dtmaxdt.strftime('%Y%b%d%H%M')) , pdf=pdf, bins=bins )
