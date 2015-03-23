from iris.coord_categorisation import add_categorised_coord
import iris
import numpy as np
import pdb

def add_hour_of_day(cube, coord, name='hour'):
    try:
        add_categorised_coord(cube, name, coord,
          lambda coord, x: coord.units.num2date(x).hour)
    except ValueError, e:
        print e
        pass

def add_day_of_year(cube, time_coords):
     try:
        iris.coord_categorisation.add_day_of_year(cube, time_coords, name='day_of_year')
     except ValueError, e:
        print e
        pass

def random_string(str_len):
    import string
    import random

    return ''.join([random.choice(string.letters) for i in xrange(str_len)])

temp_save_dir = '/nfs/a90/eepdw/temp_files/'
    
import iris
import numpy as np

# def auxcoord_flatten(cube, auxcoord):
#     """
#     Reshape cube so that auxcoord becomes 1-dimensional.
   
#     Every dimension of cube must have a DimCoord
#     """
   
#     auxcoord = cube.coord(auxcoord)

#     pdb.set_trace()
       
#     if len(auxcoord.shape) == 1:
#         return cube
   
#     auxcoord_dims = cube.coord_dims(auxcoord)
#     auxcoord_largest_dim = np.array(auxcoord.shape).argmax()
#     auxcoord_max_len = auxcoord.shape[auxcoord_largest_dim]

#     auxcoord_dim_coord_longest = cube.dim_coords[auxcoord_dims[auxcoord_largest_dim]].name()
#     auxcoord_shared_dim_coords = []
#     for dim in auxcoord_dims:
#         if dim != auxcoord_largest_dim:
#             auxcoord_shared_dim_coords += [coord.name() for coord in cube.coords(dimensions=dim)]
   
#     other_dim_coords = [cube.dim_coords[i].name() for i in range(cube.ndim) if i not in cube.coord_dims(auxcoord)]
   
#     temp_cubelist = iris.cube.CubeList()
   
#     for i, cslice in enumerate(cube.slices([auxcoord_dim_coord_longest] + other_dim_coords)):
#         demote_coord = cslice.coord(auxcoord_dim_coord_longest)
#         demote_coord_dim = cslice.coord_dims(demote_coord)[0]
#         cslice.remove_coord(demote_coord)
#         cslice.add_aux_coord(demote_coord, demote_coord_dim)
#         cslice.add_dim_coord(iris.coords.DimCoord(np.arange(auxcoord_max_len)+auxcoord_max_len*i,
#                                                   long_name='pseu
 #                                                  do'),
#                              demote_coord_dim)
       
#         for coordname in auxcoord_shared_dim_coords:
#             coord = cslice.coord(coordname)
#             cslice.remove_coord(coord)
#             points = np.repeat(coord.points, auxcoord_max_len)
#             if coord.has_bounds():
#                 bounds = np.repeat(coord.bounds, auxcoord_max_len, 0)
#             else:
#                 bounds = None
#             coord = coord_from_coord(points, coord, bounds, False)
#             cslice.add_aux_coord(coord, demote_coord_dim)
       
#         temp_cubelist.append(cslice)
         
#     return temp_cubelist.concatenate_cube()


# def coord_from_coord(points, samplecoord, bounds=None, dimcoord=True):
#     """
#     Return a coord with all metadata from samplecoord
   
#     INPUTS
   
#     points       numpy array
#                  as for iris.coord.DimCoord and iris.coord.AuxCoord
                 
#     samplecoord  iris.coord.DimCoord or iris.coord.AuxCoord
   
#     bounds       numpy array
#                  as for iris.coord.DimCoord and iris.coord.AuxCoord
                 
#     dimcoord     boolean
#                  set True to return a DimCoord, False to return an AuxCoord
#     """
#     if dimcoord:
#         coordtype = iris.coords.DimCoord
#     else:
#         coordtype = iris.coords.AuxCoord
   
#     return coordtype(points,
                     # standard_name=samplecoord.standard_name,
                     # long_name=samplecoord.long_name,
                     # var_name=samplecoord.var_name,
                     # units=samplecoord.units,
                     # bounds=bounds,
                     # attributes=samplecoord.attributes,
                     # coord_system=samplecoord.coord_system)


def FlattenIfTimeMultipleCoords(cube, exp_id):
    
    auxcoord='time'

    import os
   
    if len(cube.coord(auxcoord).points.shape)>1:
    
       import os

       other_dim_coords = [cube.dim_coords[i].name() for i in range(cube.ndim) if i not in cube.coord_dims(auxcoord)]

       #pdb.set_trace()

       save_file = '%s%s_%s_flatten.pp' \
                         % (temp_save_dir, cube.standard_name, exp_id)

       try:
           os.remove(save_file)
       except Exception, e:
           print e

       #cube_list=[]
        
       for t, time_cube in enumerate(cube.slices(other_dim_coords)):
       #for t, time_cube in enumerate(cube.slices('grid_latitude', 'grid_longitude')):

            time_cube.remove_coord('forecast_period')
            time_cube.remove_coord('forecast_reference_time')
        
            #iris.save(time_cube, 
             #          temp_save_file, append=True)

            iris.save(time_cube, save_file, append=True)

            print t
        
        #cube = iris.cube.CubeList(cube_list).merge_cube()

    #pdb.set_trace()

       #cube = iris.load(temp_save_file).merge_cube()

       #iris.save(cube, '%s%s_%s_flatten.pp' \
       #                  % (temp_save_dir, cube.standard_name, exp_id) 
        #        )

       #os.remove(temp_save_file)
    cube =iris.load_cube(save_file)
    
    return cube

def FlattenIfTimeMultipleCoordsSave(cube, exp_id):
    
    auxcoord='time'

    import os
   
    if len(cube.coord(auxcoord).points.shape)>1:
    
       import os

       other_dim_coords = [cube.dim_coords[i].name() for i in range(cube.ndim) if i not in cube.coord_dims(auxcoord)]

       #pdb.set_trace()

       temp_save_file = '%s%s%s.pp' \
                         % (temp_save_dir, cube.standard_name, random_string(8))

       #cube_list=[]
        
       for t, time_cube in enumerate(cube.slices(other_dim_coords)):
       #for t, time_cube in enumerate(cube.slices('grid_latitude', 'grid_longitude')):

            time_cube.remove_coord('forecast_period')
            time_cube.remove_coord('forecast_reference_time')
        
            #iris.save(time_cube, 
             #          temp_save_file, append=True)

            iris.save(time_cube, '%s%s_%s_flatten.pp' \
                         % (temp_save_dir, cube.standard_name, exp_id), append=True)

            print t
        
        #cube = iris.cube.CubeList(cube_list).merge_cube()

       pdb.set_trace()

       #cube = iris.load(temp_save_file).merge_cube()

       #iris.save(cube, '%s%s_%s_flatten.pp' \
       #                  % (temp_save_dir, cube.standard_name, exp_id) 
        #        )

       #os.remove(temp_save_file)
    #cube =iris.load_cube( '%s%s_%s_flatten.pp' \
                         #% (temp_save_dir, cube.standard_name, exp_id))
    
    #return cube

def CubeMeanHourlyGridLatLon(cube, experiment_id):  
    '''
 
    '''

    cube = FlattenIfTimeMultipleCoords(cube, experiment_id)

    time_coords = cube.coord('time')
    
    add_hour_of_day(cube, time_coords)
        
    means_hourly = cube.aggregated_by('hour' , iris.analysis.MEAN)

    means_cube = means_hourly.collapsed(('grid_latitude', 'grid_longitude'), iris.analysis.MEAN)  

    hours_in_data = cube.coord('hour').points

    return hours_in_data, means_cube


def CubeMeanDailyGridLatLon(cube, experiment_id):  
    '''

    '''
    

    cube = FlattenIfTimeMultipleCoords(cube, experiment_id)
        
    time_coords = cube.coord('time')
    
    add_day_of_year(cube, time_coords)

    means_daily = cube.aggregated_by('hour' , iris.analysis.MEAN)

    means_cube = means_daily.collapsed(('grid_latitude', 'grid_longitude'), iris.analysis.MEAN) 

    hours_in_data = cube.coord('hour').points

    return hours_in_data, means_cube

def CubeMeanHourly(cube, experiment_id):  
    '''
  
    '''
  
    cube = FlattenIfTimeMultipleCoords(cube, experiment_id)

    time_coords = cube.coord('time')
  
    add_hour_of_day(cube, time_coords)
        
    means_hourly = cube.aggregated_by('hour' , iris.analysis.MEAN)

    hours_in_data = cube.coord('hour').points

    return hours_in_data, means_hourly


def CubeMeanDaily(cube, experiment_id):  
    '''
 
    '''

    cube = FlattenIfTimeMultipleCoords(cube, experiment_id)
    
    time_coords = cube.coord('time')
     
    add_day_of_year(cube, time_coords)

    means_daily = cube.aggregated_by('day_of_year' , iris.analysis.MEAN)

    days_in_data = cube.coord('day_of_year').points

    return days_in_data, means_daily


def LandSeaHourlyMean(cube, lsm, experiment_id):
    '''
    Takes pp cube, pp land/sea mask cube, and cube name (only for save filename)
    and returns iris cube of hourly mean
    '''
        
    for s in ([0, 1]):
    
        #import pdb
        #pdb.set_trace()

        #temp_save_file = '%s%s%s_%s.pp' \
            #             % (temp_save_dir, cube.standard_name, random_string(8), s)

        #mean_list=iris.cube.CubeList([])
        
        

        # for t, time_cube in enumerate(cube.slices(['grid_latitude','grid_longitude'])):

        #     sm = np.where(lsm_boolean, time_cube.data, np.nan)
        #     maskedcube = np.ma.masked_array(sm, sm==np.nan)
        #     #mean_list.append(time_cube.copy(data=maskedcube))
        #     iris.save(time_cube.copy(data=maskedcube), 
        #               temp_save_file, append=True)

        #del sm, maskedcube

        #pdb.set_trace()
 
        #cube = iris.load(temp_save_file).merge_cube()

        #mean = CubeMeanHourly(cube)
        #mean2=mean.collapsed(('grid_latitude', 'grid_longitude'), iris.analysis.MEAN)
        
        lsm_boolean = (lsm.data==s)
        cube.mask=lsm_boolean

        if s==0:
            hours, mean_0 = CubeMeanHourlyGridLatLon(cube, experiment_id)
        if s==1:
            hours, mean_1 = CubeMeanHourlyGridLatLon(cube, experiment_id)
            
    return mean_0, mean_1
            
def LandSeaDailyMean(cube, lsm, experiment_id):
    '''
    Takes pp cube, pp land/sea mask cube, and cube name (only for save filename)
    and returns iris cube of daily mean
    '''

    
    for s in ([0, 1]):
    
        # sm = np.where(lsm.data==s, cube.data, np.nan)
        # maskedcube = np.ma.masked_array(sm, sm==np.nan)
        # masked_iris_cube = cube.copy(data=maskedcube)

        # mean = CubeMeanDaily(cube, experiment_id)
        # mean2=mean.collapsed(('grid_latitude', 'grid_longitude'), iris.analysis.MEAN)
        
        lsm_boolean = (lsm.data==s)
        cube.mask=lsm_boolean

        if s==0:
            hours, mean_0 = CubeMeanDailyGridLatLon(cube, experiment_id)
        if s==1:
            hours, mean_1 = CubeMeanDailyGridLatLon(cube, experiment_id)
            
    return mean_0, mean_1


def PointMeanHourly(cube):  
    '''
    Takes a cube (one grid_lat, one grid_lon) as input and 
    iterates through hourly aggregates getting mean
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
    
        try:
            means.append(cube[time_idx_for_mean[1], 
                                   time_idx_for_mean[0]].collapsed('time', 
                                        iris.analysis.MEAN).data.data[0][0])    
        except Exception:
            means.append(cube[time_idx_for_mean].collapsed('time', 
                                        iris.analysis.MEAN).data.data[0][0])    
    return hours_in_data, means

