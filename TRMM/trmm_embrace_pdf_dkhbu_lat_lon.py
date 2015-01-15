import cPickle as pickle
import numpy as np
from collections import defaultdict
from netCDF4 import Dataset
from scipy.interpolate import griddata

from shapely.geometry import Point, Polygon



# Min and max lats lons from smallest model domain (dkbhu) - see spreadsheet

# Unrotated - can only be used (without unrotating first) for djznw  - global model

latmin=-6.79
latmax=33.04
lonmin=64.12
lonmax=101.87

pcp_dom, longitude_dom, latitude_dom, time_dom, time_hour = pickle.load(open('/nfs/a90/eepdw/Data/Saved_data/TRMM/trmm_emb_time_update_large.p', 'rb'))

# Load land sea mask.  TRMM land sea mask is in % of water coverage so 100% is all water

nc = Dataset('/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/TMPA_mask.nc')

#  Regrid lsm to data grid (offset b 0.125 degrees

lsm_lons, lsm_lats = np.meshgrid(nc.variables['lon'][:],nc.variables['lat'][:])
lons_data, lats_data = np.meshgrid(longitude_dom[0], latitude_dom[0])
lsm_regrid = griddata((lsm_lats.flatten(), lsm_lons.flatten()), nc.variables['landseamask'][:].flatten(), (lats_data,lons_data), method='linear')

#  Find points that are within defined polygon

lat_max_idx = np.max(np.where((latmin<=latitude_dom[0]) & (latitude_dom[0]<=latmax)))
lat_min_idx = np.min(np.where((latmin<=latitude_dom[0]) & (latitude_dom[0]<=latmax)))

lon_max_idx = np.max(np.where((lonmin<=longitude_dom[0]) & (longitude_dom[0]<=lonmax)))
lon_min_idx = np.min(np.where((lonmin<=longitude_dom[0]) & (longitude_dom[0]<=lonmax)))

pcp_dom_2 = pcp_dom[:, lat_min_idx:lat_max_idx, lon_min_idx:lon_max_idx]
lsm = lsm_regrid[lat_min_idx:lat_max_idx, lon_min_idx:lon_max_idx]

print pcp_dom.shape
print pcp_dom_2.shape
print lsm.shape
lsm=np.resize(lsm, pcp_dom_2.shape)[0]

####################################################

bins=np.linspace(0.,200., 200)

for ls in ['land', 'sea']:
    if ls=='land':
        pdf=[]
        bmin=0.
        for b in bins:
            pdf.append(((np.resize(np.abs(100-lsm), pcp_dom_2.shape)[(pcp_dom_2 <= b) & (pcp_dom_2>=bmin)])/100).sum())
            if bmin==0.:
                bmin=0.0000000001
            else:
                bmin=b

    elif ls=='sea':
        pdf=[]
        bmin=0.
        for b in bins:
            pdf.append(((np.resize(lsm, pcp_dom_2.shape)[(pcp_dom_2 <= b) & (pcp_dom_2>=bmin)])/100).sum())
            if bmin==0.:
                bmin=0.0000000001
            else:
                bmin=b
#####################################################


# Save


    np.savez("/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/Rainfall_Intensity/EMBRACERainfall_Intensity_PDF_dkbhu_latlon_%s" % ls , pdf=pdf, bins=bins )
