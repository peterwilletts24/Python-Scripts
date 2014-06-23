from iris.coords import DimCoord
from iris.cube import Cube
import iris
import numpy as np

def simple(pressure):
    c = Cube(np.zeros((20, 30), 'f4'))
    c.add_dim_coord(DimCoord(range(20), 'latitude', units='degrees'), 0)
    c.add_dim_coord(DimCoord(range(30), 'longitude', units='degrees'), 1)
    c.add_aux_coord(DimCoord(pressure, long_name='pressure', units='hPa'))
    return c

if __name__ == '__main__':
    print simple(100)

    iris.save(simple(100), 'test.pp', append=True)
    iris.save(simple(150), 'test.pp', append=True)
    for cube in iris.load('test.pp'):
        print cube
        print cube.coord('pressure')
