experiment_ids = ['dklyu']
#experiment_ids = ['dklyu', 'dkmbq', 'djzny', 'djznq', 'djzns', 'dkjxq', 'dklwu', 'dklzq', 'dkbhu' ] # All 12
experiment_ids = ['dkbhu' ]

station_list=['Calcutta / Dum Dum', 'Port Blair', 'Nagpur Sonegaon']
station_list=['Port Blair', 'Nagpur Sonegaon']
#station_list=['Patna', 'Port Blair', 'Nagpur Sonegaon', 'Gorakhpur', 'Allahabad_Bamhrauli'] 
station_list=['Patna', 'Gorakhpur', 'Allahabad_Bamhrauli']
import iris
from iris.coord_categorisation import add_categorised_coord

def add_hour_of_day(cube, coord, name='hour'):
    add_categorised_coord(cube, name, coord,
          lambda coord, x: coord.units.num2date(x).hour)
    
def PointMeanHourly(cube):  
    '''
    Takes a cube as input and iterates through hourly aggregates getting mean
    Returns a list of means and another of corresponding hours
    Written as a work around for when time is in 2 coordinates and iris functions
    will not work
    '''
    
    time_coords = cube.coord('time')
    try:
        add_hour_of_day(cube, time_coords)
    except ValueError, e:
        print e
        pass

    hours_in_data = np.unique(cube.coord('hour').points)
        
    means=[]
    
    for h in hours_in_data:
        #print h
    
        time_idx_for_mean = np.where(cube.coord('hour').points==h)
            #print mslp_point[time_idx_for_mean]
    
        means.append(mslp_point[time_idx_for_mean[1], 
                                   time_idx_for_mean[0]].collapsed('time', 
                                        iris.analysis.MEAN).data.data[0][0])    
    return hours_in_data, means

import numpy as np
import datetime

wmo_or_icao = 'wmo'

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

        mslp = iris.load_cube('/nfs/a90/eepdw/Data/EMBRACE/%s/%s/16222.pp' % (expmin1, experiment_id))
        
        cs = mslp.coord_system('CoordSystem')
                    
        #lons, lats = np.meshgrid(lon, lat) 
        grid_lon, grid_lat = iris.analysis.cartography.rotate_pole\
                                (stat_lon,stat_lat, cs.grid_north_pole_longitude, cs.grid_north_pole_latitude)
            
        grid_lon = grid_lon % 360
        
        sample_points = [('grid_latitude', grid_lat),
                         ('grid_longitude', grid_lon)]
        
        mslp_point = iris.analysis.interpolate.linear(mslp, sample_points)
        
        hours, means = PointMeanHourly(mslp_point)
        
        np.savez('/nfs/a90/eepdw/Data/EMBRACE/%s/%s/%s_same_lat_lon_as_%s_hourly_mean_mslp' 
                % (expmin1, experiment_id, experiment_id, stat.replace(' / ', '_').replace(' ', '_')), hours=hours, means=means )
        
