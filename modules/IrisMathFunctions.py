import iris
import numpy as np

def SpHumToMixingRatio(cube):
    '''
    Take iris cube of specific humidity as input and converts to mixing ratio
    Rearrange this equation
    http://glossary.ametsoc.org/wiki/Specific_humidity
    '''

    denominator = iris.analysis.maths.subtract(cube, 1.)*-1.

    return cube/denominator

def TotalColumnWaterVapour(mixing_ratio, pressure):
    '''
    Take iris cubes of mixing ratio and pressure (408 in standard UM diagnostics)
    and return total column water vapour cube
    '''

    density_of_water = 1./100. # 1 g/cm^3

    for c, coord in enumerate(mixing_ratio.coords()):
        if coord.standard_name=='model_level_number':
            mod_lev_num_idx = c
    
    tcwv = (1/(density_of_water*9.81))*np.trapz(mixing_ratio.data, pressure.data, axis=mod_lev_num_idx)
    tcwv = mixing_ratio[:,:,0,:,:].copy(data=tcwv)

    return tcwv

    
