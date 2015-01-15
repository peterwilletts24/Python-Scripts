import cPickle as pickle
import numpy as np
from collections import defaultdict
from netCDF4 import Dataset
from scipy.interpolate import griddata

from shapely.geometry import Point, Polygon


polygon = Polygon(((73., 21.), (83., 16.), (87., 22.), (75., 27.)))

#lon_high = 116 
#lon_low = 30.5

#lat_high= 40
#lat_low=-11.25

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

bins=np.linspace(0,200, 200)
pdf=[]
bmin=0.
for b in bins:
    pdf.append(pcp_dom_2[(pcp_dom_2 <= b) & (pcp_dom_2>=bmin)].shape[0])
    if bmin==0:
        bmin=0.0000000001
    else:
        bmin=b
#####################################################


# Save


np.savez("/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/Rainfall_Intensity/EMBRACE_Monsoon_Trough_Rainfall_Intensity_PDF" , pdf=pdf, bins=bins )
