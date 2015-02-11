#experiment_ids = ['dklyu']
experiment_ids = ['djzny', 'dklyu', 'dkmbq', 'djznq', 'djzns', 'dkjxq', 'dklwu', 'dklzq'] # All 12
#experiment_ids = ['dklwu', 'dklzq', 'dkbhu' ]

#station_list=['Calcutta / Dum Dum', 'Port Blair', 'Nagpur Sonegaon', 'Gorakhpur', 'Allahabad / Bamhrauli', 'Patna']
#station_list=['Varanasi / Babatpur']
station_list=['Bombay / Santacruz']
# cube_names=['x_wind','y_wind']
# cube_names_2=['eastward_wind', 'northward_wind']
# No standard or long names for 10m wind

import iris
import numpy as np
import datetime

import imp
imp.load_source('IrisTimeMeanRoutines', '/nfs/see-fs-01_users/eepdw/python_scripts/modules/IrisTimeMeanFunctions.py')
from IrisTimeMeanRoutines import *

wmo_or_icao = 'wmo'

header = np.loadtxt('/nfs/a90/eepdw/Data/Midas/GL_Column_Headers.csv', dtype=str, delimiter=',')
header = [h.strip(' ') for h in header]

plot_list = np.load('/nfs/a90/eepdw/Data/Midas/station_plot_list_%s.npy' % wmo_or_icao)
plot_list=np.array([['Bombay / Santacruz',  '19.116666666666667', '72.85', '43003']])
#p_level=925.

for stat in station_list:
    pdb.set_trace()
    station_info = plot_list[np.where(plot_list==stat)[0]][0]
    
    stat_lat = np.array(float(station_info[1]))
    stat_lon = np.array(float(station_info[2]))
    
    print stat
    
    for experiment_id in experiment_ids:
        
        print experiment_id
    
        expmin1 = experiment_id[:-1]
        
        wind_file = '/nfs/a90/eepdw/Data/EMBRACE/%s/%s/3209.pp' % (expmin1, experiment_id)
    
        #try:
    
        east_wind = iris.load(wind_file)[0]
        north_wind = iris.load(wind_file)[1]      
                            
        #except iris.exceptions.ConstraintMismatchError: 

            #east_wind = iris.load_cube(wind_file, cube_names_2[0] & iris.Constraint(pressure=p_level))
            #north_wind = iris.load_cube(wind_file, cube_names_2[1] & iris.Constraint(pressure=p_level))      
                                
        cs = east_wind.coord_system('CoordSystem')
                    
        #lons, lats = np.meshgrid(lon, lat) 
        grid_lon, grid_lat = iris.analysis.cartography.rotate_pole\
                                (stat_lon,stat_lat, cs.grid_north_pole_longitude, cs.grid_north_pole_latitude)
            
        grid_lon = grid_lon % 360
        
        sample_points = [('grid_latitude', grid_lat),
                         ('grid_longitude', grid_lon)]
        
        east_wind_point = iris.analysis.interpolate.linear(east_wind, sample_points)
        north_wind_point = iris.analysis.interpolate.linear(north_wind, sample_points)
        
        #pdb.set_trace()

        east_hours, east_means = PointMeanHourly(east_wind_point)
        north_hours, north_means = PointMeanHourly(north_wind_point)
        
        #pdb.set_trace()

        np.savez('/nfs/a90/eepdw/Data/EMBRACE/%s/%s/%s_same_lat_lon_as_%s_hourly_mean_wind_10m' \
                % (expmin1, experiment_id, experiment_id, stat.replace(' / ', '_').replace(' ', '_')), \
                 east_hours=east_hours, east_means=east_means, 
                 north_hours=north_hours, north_means=north_means)
