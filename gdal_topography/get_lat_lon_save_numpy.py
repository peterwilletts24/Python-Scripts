from osgeo import osr, gdal
import numpy as np

file="/nfs/a90/eepdw/Data/Observations/Satellite/GTOPO3_1km_topography/gt30e060n40.tif"
# get the existing coordinate system
ds = gdal.Open('%s' % file)
#cs= osr.SpatialReference()
#cs.ImportFromWkt(ds.GetProjectionRef())
 
#get the point to transform, pixel (0,0) in this case
width = ds.RasterXSize
height = ds.RasterYSize
gt = ds.GetGeoTransform()
minx = gt[0]
miny = gt[3] + width*gt[4] + height*gt[5] 
maxx = gt[0] + width*gt[1] + height*gt[2]
maxy = gt[3]

latPxSz = gt[5]
lonPxSz = gt[1]

lons = np.arange(minx, maxx, lonPxSz)
lats = np.arange(maxy, miny, latPxSz)

data=np.array(gdal.Open("%s" % file).ReadAsArray())

np.savez("/nfs/a90/eepdw/Data/Observations/Satellite/GTOPO3_1km_topography/gt30e060n40", lats_gtopo30=lats, lons_gtopo30=lons, gtopo30_data=data)
