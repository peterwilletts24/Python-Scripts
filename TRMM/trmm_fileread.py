import re
from collections import namedtuple
import numpy as np

trmm_files = np.zero(1, dtype = [('year',a4), 
                                  ('month',a2), 
                                  ('day',a2), 
                                  ('hour',a2)])

with open('/nfs/see-fs-01_users/eepdw/python_scripts/filenamelist/trmm_netcdffilelist', 'r') as f:
    netcdf_filelist=f.readlines()

ncdfl2 = [elem.strip().split("/") for elem in netcdf_filelist]

#
#trmm_files = namedtuple("trmm_files", "year month day hour")
#print s2[1:51]
#print netcdf_filelist[1]
#print ncdfl2[0][9]
#print len(ncdfl2)

for i, item in enumerate(ncdfl2):
    trmm_files.year[0][i]=  item[9]

print trmm_files.year
#/nfs/a80/earceb/model_data/observations/satellite/TRMM/2011
# nappy.convertNAToCSV(na_file, annotation=True)
