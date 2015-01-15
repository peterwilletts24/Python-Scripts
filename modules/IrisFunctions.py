import iris.unit as unit
import matplotlib

import datetime

import numpy as np

import imp

"""
imp.load_source('IrisFunctions', '/nfs/see-fs-01_users/eepdw/python_scripts/modules/IrisFunctions.py')
from IrisFunctions import *
"""

def ConvertHoursSince1970ToString(hours_since, string_format):
    """ Takes numpy array of 'hours since 1970' and string format
     and returns numpy array of date strings"""

    u = unit.Unit('hours since 1970-01-01 00:00:00',calendar='gregorian')

    return np.array([u.num2date(da).strftime(string_format) for da in np.array(hours_since)])

def ConvertHoursSince1970ToDatetime(hours_since):
    """ Takes numpy array of 'hours since 1970' and string format
     and returns numpy array of date strings"""

    u = unit.Unit('hours since 1970-01-01 00:00:00',calendar='gregorian')

    return np.array([u.num2date(da) for da in np.array(hours_since)])

def ConvertTimeStampStringToDatetime(time_string_array, string_format):
    """ Takes numpy array of time stamps  in string format and format of string
     (e.g. '%Y-%m-%d %H') and returns numpy array of datetimes"""

    return np.array([datetime.datetime.strptime(da, string_format) for da in time_string_array] )
#def 


# Function to get index of cube coords

# Divergence

def Divergence(u_wind, v_wind):
    """ 
    Takes u and v wind iris cubes (on latitude and longitude grdi)
    as input and returns divergence of wind fields
    http://www.met.wau.nl/education/atmospract/unit19/Div_vortUK.pdf - Appendix A
    TO DO add assertions for same pressure level, etc
    """
    r = 6371.22e6  #Average radius of earth in metres, same as ERA Interim

    dudlon = iris.analysis.calculus.differentiate(u_wind, 'longitude')
    dudlon.units=None

    cos_lats=iris.analysis.cartography.cosine_latitude_weights(v_wind)

    v_cos_lats_interp = iris.analysis.interpolate.regrid(iris.analysis.maths.multiply(v_wind, cos_lats), dudlon, mode='bilinear')

    dvdlat = iris.analysis.calculus.differentiate(v_cos_lats_interp, 'latitude')
    dvdlat.units=None

    if dudlon.shape != dvdlat.shape:
        dvdlat.transpose([1,0,2,3])


    second_term=iris.analysis.maths.add(dudlon, dvdlat)

    first_term=1/(r*cos_lats)   

    return iris.analysis.maths.multiply(second_term, first_term) 

def Vorticity(u_wind, v_wind):
    """ 
    Takes u and v wind iris cubes (on latitude and longitude grdi)
    as input and returns vorticity of wind fields
    http://www.met.wau.nl/education/atmospract/unit19/Div_vortUK.pdf - Appendix A
    TO DO add assertions for same pressure level, etc
    """
    r = 6371.22e6  #Average radius of earth in metres, same as ERA Interim

    dvdlon = iris.analysis.calculus.differentiate(v_wind, 'longitude')
    dvdlon.units=None

    cos_lats=iris.analysis.cartography.cosine_latitude_weights(u_wind)

    u_cos_lats_interp = iris.analysis.interpolate.regrid(iris.analysis.maths.multiply(u_wind, cos_lats), dudlon, mode='bilinear')

    dudlat = iris.analysis.calculus.differentiate(u_cos_lats_interp, 'latitude')
    dudlat.units=None

    if dudlon.shape != dvdlat.shape:
        dudlat.transpose([1,0,2,3])


    second_term=iris.analysis.maths.subtract(dvdlon, dudlat)

    first_term=1/(r*cos_lats)   

    return iris.analysis.maths.multiply(second_term, first_term) 
    
    
