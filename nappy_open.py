import os
import sys
#import cdms
import glob
import nappy

rad_flist = glob.glob ('/nfs/a80/eepdw/india_radiosonde_2011/*/2011/*_20110[89].na')

for i in rad_flist:

 linestr=i.rstrip()

 f=nappy.openNAFile(linestr)
 print f
