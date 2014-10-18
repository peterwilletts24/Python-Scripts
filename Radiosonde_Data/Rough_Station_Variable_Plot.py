########################################
# Read radiosonde files from Integrated Global RadiosondeArchive
#, separate soundings and save in python
#
# http://www1.ncdc.noaa.gov/pub/data/igra
# Created by: Peter Willetts
# Created on: 25/06/2014
#
########################################
#
#
###################################################

import glob

rad_flist = glob.glob ('/nfs/a90/eepdw/Data/Observations/Radiosonde_downloaded_from_NOAA_GUAN/derived_parameters/*.dat')
