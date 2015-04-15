# Now make a simple example using the custom projection.
import pdb

import sys
import os
import pkg_resources
pkg_resources.require('matplotlib==1.4.0')
import datetime
from dateutil.relativedelta import relativedelta
import re
import math

from matplotlib.ticker import ScalarFormatter, MultipleLocator
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt
from StringIO import StringIO
import numpy as np

from numpy import load

# Exception handling, with line number and stuff
import linecache
import sys

def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print 'EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj)

import imp
imp.load_source('SoundingRoutines', '/nfs/see-fs-01_users/eepdw/python_scripts/Tephigram/Sounding_Routines.py')
imp.load_source('TephigramPlot', '/nfs/see-fs-01_users/eepdw/python_scripts/Tephigram/Tephigram_Functions.py')
from TephigramPlot import *
from SoundingRoutines import *

imp.load_source('GeogFuncs', '/nfs/see-fs-01_users/eepdw/python_scripts/modules/GeogFunctions.py')
from GeogFuncs import *

pmin=200.

station_list_cs=[42182, 43003, 43014, 42867, 43371, 43353, 43285,  43192, 43150, 42339, 40990, 40948]
#station_list_cs=[43003]
date_min=datetime.datetime(1960,5,1,0,0,0)
date_max=datetime.datetime(2014,10,1,0,0,0)

delta = relativedelta(weeks=+1)    
        
variable_list={'pressures': 0, 'temps':1, 'dewpoints':2, 'winddirs':3, 'windspeeds':4, 'pot_temp':5, 
               'sat_vap_pres':6, 'vap_press':7, 'rel_hum':8, 'wvmr':9, 'sp_hum':10, 'sat_temp':11, 
               'theta_e':12, 'theta_e_sat':13, 'theta_e_minus_theta_e_sat':14}

variable_list_line={'lcl_temp': 0, 'lcl_vpt':1, 'pbl_pressure':2, 'surface_pressure':3, 'T_eq_0':4}

def variable_name_index_match(variable, variable_list):
    for key, value in variable_list.iteritems():   # iter on both keys and values
        if key.startswith('%s' % variable) and key.endswith('%s' % variable):
            arr_index_var=value 
    return arr_index_var



    # Parse the data
for stat in station_list_cs:

    station_name,la,lo, st_height=StationInfoSearch(stat)
    
    load_file = load('/nfs/a90/eepdw/Data/Observations/Radiosonde_Numpy/Radiosonde_Cross_Section_'
                           'IND_SOUNDING_INTERP_MEAN_Climat_%s_%s_%s_%s.npz'
                           % (date_min.strftime('%Y%m%d'), date_max.strftime('%Y%m%d'), delta, stat))
    data=load_file['date_bin_mean_all_dates_one_station']

    dates=load_file['dates_for_plotting']
    for bin in range(data.shape[0]):
        try:           
            p=data[bin,0,:]/100
            T=data[bin,1,:]-273.15
            Td=T-data[bin,2,:]
            h=data[bin,15,:]
            da=dates[bin]
            #print T
            #print p
            #print Td

            #pdb.set_trace()

            #u_wind,v_wind = u_v_winds(data[bin,3,:], data[bin,4,:])
            u_wind,v_wind = data[bin,-2,:], data[bin,-1,:]
            # Create a new figure. The dimensions here give a good aspect ratio
            fig = plt.figure(figsize=(10, 8), frameon=False)
            #fig.patch.set_visible(False)
            
            tephigram_plot_height=0.85
            tephigram_plot_bottom=.085
            ax = fig.add_axes([.085,tephigram_plot_bottom,.65,tephigram_plot_height], projection='skewx', frameon=False, axisbg='w')
            ax.set_yscale('log')
            plt.grid(True)

 

            #pdb.set_trace()
            tmax=math.ceil(nanmax(T)/10)*10
            tmin=math.floor(nanmin(Td[p>400])/10)*10
            pmax=math.ceil(nanmax(p)/50)*50

            P=linspace(pmax,pmin,37)

            w = array([0.0001,0.0004,0.001, 0.002, 0.004, 0.007, 0.01, 0.016, 0.024, 0.032, 0.064, 0.128])
            ax.add_mixratio_isopleths(w,linspace(pmax, 700., 37),color='m',ls='-',alpha=.5,lw=0.5)
            ax.add_dry_adiabats(linspace(-40,40,9),P,color='k',ls='-',alpha=.5,lw=0.8)
            ax.add_moist_adiabats(linspace(-40,40,18),P,color='k',ls='--',alpha=.5,lw=0.8, do_labels=False)
            ax.other_housekeeping(pmax, pmin, 40,-40) 

            wbax = fig.add_axes([0.75,tephigram_plot_bottom,0.12,tephigram_plot_height],frameon=False, sharey=ax, label='barbs')
            ax_text_box = fig.add_axes([0.85,0.085,.12,tephigram_plot_height], frameon=False, axisbg='w')

          # Plot the data using normal plotting functions, in this case using semilogy
 
            ax.semilogy(T, p, 'r', linewidth=2)
            ax.semilogy(Td, p, 'r',linewidth=2)

            # row_labels=(
            #             'SLAT',
            #             'SLON',
            #             'SELV',
            #             'SHOW',
            #             'LIFT',
            #             'LFTV',
            #             'SWET',
            #             'KINX',
            #             'CTOT',
            #             'VTOT',
            #             'TOTL',
            #             'CAPE',
            #             'CINS',
            #             'CAPV',
            #             'CINV',
            #             'LFCT',
            #             'LFCV',
            #             'BRCH',
            #             'BRCV',
            #             'LCLT',
            #             'LCLP',
            #             'MLTH', 
            #             'MLMR',
            #             'THCK',
            #             'PWAT')

            # variable='pbl_pressure'
            # var_index = variable_name_index_match(variable, variable_list_line)
            # print load_file['date_bin_mean_all_dates_one_station_single'].shape
            # pbl_pressure =  load_file['date_bin_mean_all_dates_one_station_single'][bin,0,var_index]

            # print pbl_pressure

            # EQLV, pp, lclp,lfcp, lclt, delta_z, CAPE, CIN=CapeCinPBLInput(p, T, Td, h, st_height, pbl_pressure/100)

            # print lclp


            # table_vals=(
            #     #'%s' % station_name,
            #     #'Climatology - Week beg. %s' % da, 
            #     '%s' % la,
            #     '%s' % lo,
            #     '%s' % st_height,
            #     '%.1f' % ShowalterIndex(T, Td, p),                       #      ['Showalter index',
            #     '%.1f' % LiftedIndex(T, Td, p, h, st_height),          #     'Lifted index', 
            #     '--',                                                                    #    'LIFT computed using virtual temperature',
            #     '--',                                                                    #   'SWEAT index',
            #     '%.1f' % KIndex(T, Td, p),          #    'K index',
            #     '%.1f' % CrossTotalsIndex(T, Td, p),                   #   'Cross totals index',
            #     '%.1f' % VerticalTotalsIndex(T, p),                      #   'Vertical totals index',
            #     '%.1f' % TotalTotalsIndex(T, Td, p),                     #   'Total totals index',
            #     '%.1f' % CAPE,                                                   #   'CAPE',
            #     '%.1f' % CIN,                                                      #   'CIN',
            #     '--',                                                                    #   'CAPE using virtual temperature',
            #     '--',                                                                    #   'CINS using virtual temperature',
            #     '%.1f' % lfcp,                                                     #   'Level of free convection',
            #     '--',                                                                    #   'LFCT using virtual temperature',
            #     '--' ,                                                                   #   'Bulk Richardson number',
            #     '--',                                                                    #   'Bulk richardson using CAPV',
            #     '%.1f' % lclt,                                                        #  'Temp [K] of the Lifted Condensation Level',
            #     '%.1f' % lclp ,                                                       #  'Pres [hPa] of the Lifted Condensation Level',
            #     '--',                                                                    #  'Mean mixed layer potential temperature', 
            #     '--',                                                                    #  'Mean mixed layer mixing ratio',
            #     '--',                                                                    #  '1000 hPa to 500 hPa thickness',
            #     '--')                                                                    #  'Precipitable water [mm] for entire sounding'] 
            
            # Wind barbs

            barbs_idx=np.logspace(np.log10(10),np.log10(max(len(u_wind))),num=32).astype(int)

            wbax.set_yscale('log')
            wbax.xaxis.set_ticks([],[])
            wbax.yaxis.grid(True,ls='-',color='y',lw=0.5)

            wbax.set_xlim(-1.5,1.5)
            wbax.get_yaxis().set_visible(False)
            wbax.set_ylim(pmax+100,pmin)
                
            wbax.barbs((zeros(p.shape))[barbs_idx-1],p[barbs_idx-1], u_wind[barbs_idx-1], v_wind[barbs_idx-1])

            # Disables the log-formatting that comes with semilogy
            ax.yaxis.set_major_formatter(ScalarFormatter())
            ax.set_yticks(linspace(100,1000,10))
            ax.set_ylim(pmax,pmin)

            ax.set_xlim(-40.,40.)
            ax.xaxis.set_ticks([],[])
     
            ax_text_box.xaxis.set_visible(False)
            ax_text_box.yaxis.set_visible(False)
            for tick in wbax.yaxis.get_major_ticks():
                # tick.label1On = False
                pass
                #wbax.get_yaxis().set_tick_params(size=0,color='y')
                

            # y_loc=1.

            # max_string_length = max([len(line) for line in row_labels])
        
            # for t,r in zip(row_labels,table_vals):
            #     label_rightjust=('{:>%i}' % max_string_length).format(t)
      
            #     ax_text_box.text(0.5, y_loc, ' %s:' % (label_rightjust), size=8, horizontalalignment='right')
            #     ax_text_box.text(0.5, y_loc, ' %s' % (r), size=8, horizontalalignment='left')
            #     y_loc-=0.04

            fig.text(.02,0.965, '%s  %s' %(stat, station_name), size=12, horizontalalignment='left')
            fig.text(.02,0.035, 'Climatology - Week beg. %s ' %(da.strftime('%m-%d')), size=12, horizontalalignment='left')
            #plt.show()
            plt.savefig('/nfs/a90/eepdw/Figures/Radiosonde/Tephigrams/Weekly_Climatology/Weekly_Climatology_%s_%s_%s_Skew_T.png' % (station_name.replace('/','_').replace(' ', '_'), stat, da.strftime('%Y%m%d')))
            plt.close()
        except Exception:
            print PrintException()
