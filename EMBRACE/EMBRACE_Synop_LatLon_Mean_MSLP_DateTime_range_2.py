experiment_ids = ['dklyu']
experiment_ids = ['dklyu', 'dkmbq', 'djzny', 'djznq', 'djzns', 'dkjxq', 'dklwu', 'dklzq', 'dkbhu' ] # All 12

station_list=['Calcutta / Dum Dum', 'Port Blair', 'Nagpur Sonegaon']
#station_list=['Port Blair', 'Nagpur Sonegaon']
#station_list=['Calcutta / Dum Dum', 'Gorakhpur']

station_list=['Patna', 'Port Blair', 'Nagpur Sonegaon', 'Gorakhpur', 'Allahabad_Bamhrauli', 'Calcutta / Dum Dum' ] 

#station_list=['Patna', 'Gorakhpur', 'Allahabad_Bamhrauli']

import iris
import iris.unit as unit
from iris.coord_categorisation import add_categorised_coord
import numpy as np
import datetime

import pdb

import imp
imp.load_source('IrisTimeMeanRoutines', '/nfs/see-fs-01_users/eepdw/python_scripts/modules/IrisTimeMeanFunctions.py')
from IrisTimeMeanRoutines import *

wmo_or_icao = 'wmo'

dtmindt = datetime.datetime(2011,8,19,0,0,0)
dtmaxdt = datetime.datetime(2011,8,24,0,0,0)

dtmin = unit.date2num(dtmindt, 'hours since 1970-01-01 00:00:00', unit.CALENDAR_STANDARD)
dtmax = unit.date2num(dtmaxdt, 'hours since 1970-01-01 00:00:00', unit.CALENDAR_STANDARD)
time_constraint = iris.Constraint(time= lambda t: dtmin <= t.point <= dtmax)

header = np.loadtxt('/nfs/a90/eepdw/Data/Midas/GL_Column_Headers.csv', dtype=str, delimiter=',')
header = [h.strip(' ') for h in header]

plot_list = np.load('/nfs/a90/eepdw/Data/Midas/station_plot_list_%s.npy' % wmo_or_icao)

for stat in station_list:
    station_info = plot_list[np.where(plot_list==stat)[0]][0]
    
    stat_lat = np.array(float(station_info[1]))
    stat_lon = np.array(float(station_info[2]))
    
    print stat
    
    for experiment_id in experiment_ids:
        
        print experiment_id
    
        expmin1 = experiment_id[:-1]

        mslp = iris.load_cube('/nfs/a90/eepdw/Data/EMBRACE/%s/%s/16222.pp' % (expmin1, experiment_id), time_constraint)
        
        cs = mslp.coord_system('CoordSystem')
                    
        #lons, lats = np.meshgrid(lon, lat) 
        grid_lon, grid_lat = iris.analysis.cartography.rotate_pole\
                                (stat_lon,stat_lat, cs.grid_north_pole_longitude, cs.grid_north_pole_latitude)
            
        grid_lon = grid_lon % 360
        
        sample_points = [('grid_latitude', grid_lat),
                         ('grid_longitude', grid_lon)]
        
        mslp_point = iris.analysis.interpolate.linear(mslp, sample_points)
        
        hours, means = PointMeanHourly(mslp_point)
        
        np.savez('/nfs/a90/eepdw/Data/EMBRACE/%s/%s/%s_same_lat_lon_as_%s_hourly_mean_mslp_%s_%s' 
                 % (expmin1, experiment_id, experiment_id, stat.replace(' / ', '_').replace(' ', '_'), 
                    dtmindt.strftime('%Y%b%d%H%M'), dtmaxdt.strftime('%Y%b%d%H%M'))
                   , hours=hours, means=means)
        
