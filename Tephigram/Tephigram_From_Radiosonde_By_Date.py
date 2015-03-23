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

import imp
imp.load_source('SoundingRoutines', '/nfs/see-fs-01_users/eepdw/python_scripts/Tephigram/Sounding_Routines.py')
imp.load_source('TephigramPlot', '/nfs/see-fs-01_users/eepdw/python_scripts/Tephigram/Tephigram_Functions.py')
from TephigramPlot import *
from SoundingRoutines import *

imp.load_source('GeogFuncs', '/nfs/see-fs-01_users/eepdw/python_scripts/modules/GeogFunctions.py')
from GeogFuncs import *

pmin=200.

station_list_cs=[43003, 43014, 42867, 43371, 43353, 43285,  43192, 43150, 42339, 40990, 40948]
#station_list_cs=[42867]
#station_list_cs=[42971, 42339]
#station_list_cs=[42809]

plot_dates = [datetime.datetime(2012,5,16,0,0,0), 
                     datetime.datetime(2012,6,1,0,0,0), 
                     datetime.datetime(2012,6,15,0,0,0), 
                     datetime.datetime(2012,7,15,0,0,0)]


#plot_dates = [datetime.datetime(2012,5,15,0,0,0),
         #            datetime.datetime(2012,5,16,0,0,0), 
          #           datetime.datetime(2012,5,17,0,0,0), 
            #         datetime.datetime(2012,5,18,0,0,0), 
           #          datetime.datetime(2012,5,19,0,0,0)]

plot_dates = [datetime.datetime(2012,5,31,0,0,0),
                    datetime.datetime(2012,6,1,0,0,0), 
                     datetime.datetime(2012,6,2,0,0,0), 
                     datetime.datetime(2012,6,3,0,0,0), 
                    datetime.datetime(2012,6,5,0,0,0),
                    datetime.datetime(2012,7,1,0,0,0)]

#plot_dates = [datetime.datetime(2012,6,14,0,0,0),
#                     datetime.datetime(2012,6,15,0,0,0), 
 #                    datetime.datetime(2012,6,16,0,0,0), 
  #                   datetime.datetime(2012,6,17,0,0,0), 
   #                  datetime.datetime(2012,6,18,0,0,0)]

#plot_dates = [datetime.datetime(2012,7,14,0,0,0),
   #                  datetime.datetime(2012,7,15,0,0,0), 
   #                  datetime.datetime(2012,7,16,0,0,0), 
     #                datetime.datetime(2012,7,17,0,0,0), 
     #                datetime.datetime(2012,7,18,0,0,0)]

plot_dates = [datetime.datetime(2012,6,2,0,0,0), 
                     datetime.datetime(2012,6,3,0,0,0), 
                    datetime.datetime(2012,6,5,0,0,0),
                    datetime.datetime(2012,7,1,0,0,0)]

plot_dates = [datetime.datetime(2012,5,15,0,0,0),
                     datetime.datetime(2012,5,16,0,0,0), 
                     datetime.datetime(2012,5,17,0,0,0), 
                    datetime.datetime(2012,5,18,0,0,0), 
                    datetime.datetime(2012,5,19,0,0,0),
                    datetime.datetime(2012,5,31,0,0,0),
                    datetime.datetime(2012,6,1,0,0,0), 
                     datetime.datetime(2012,6,2,0,0,0), 
                     datetime.datetime(2012,6,3,0,0,0), 
                    datetime.datetime(2012,6,5,0,0,0),
                    datetime.datetime(2012,6,14,0,0,0),
                    datetime.datetime(2012,6,15,0,0,0), 
                     datetime.datetime(2012,6,16,0,0,0), 
                     datetime.datetime(2012,6,17,0,0,0), 
                     datetime.datetime(2012,6,18,0,0,0),
                    datetime.datetime(2012,7,14,0,0,0),
                    datetime.datetime(2012,7,15,0,0,0), 
                     datetime.datetime(2012,7,16,0,0,0), 
                     datetime.datetime(2012,7,17,0,0,0), 
                     datetime.datetime(2012,7,18,0,0,0),
                     datetime.datetime(2012,6,2,0,0,0), 
                     datetime.datetime(2012,6,3,0,0,0), 
                    datetime.datetime(2012,6,5,0,0,0),
                    datetime.datetime(2012,7,1,0,0,0)]


date_min=datetime.datetime(2011,5,1,0,0,0)
date_max=datetime.datetime(2014,10,1,0,0,0)

match_header = re.compile(r'(#.....20|#.....19)')

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


    #pd_year=  [p.year for p in plot_date]
    #pd_month=  [p.month for p in plot_dates]
    #pd_day=  [p.day for p in plot_dates]

    station_name,la,lo, st_height = StationInfoSearch(stat)
    
                
    load_file=load('/nfs/a90/eepdw/Data/Observations/Radiosonde_Numpy/Radiosonde_Single_Station_PRESSURES__IND_INTERP_SOUNDING_%s.npz' % stat)

    #pressures_for_plotting=pressures_for_plotting, single_value_vars=single_value_vars, dates_for_plotting_single=dates_for_plotting_single, dates_for_plotting_single_single=dates_for_plotting_single_single, variable_list=variable_list, variable_list_line=variable_list_line)
    #pdb.set_trace()

    data=load_file['pressures_for_plotting']

    dates=load_file['dates_for_plotting_single']

    d_year = [d.year for d in dates]
    d_month = [d.month for d in dates]
    d_day = [d.day for d in dates]

    dates_single = load_file['dates_for_plotting_single_single']

    d_year_single = [d.year for d in dates_single]
    d_month_single = [d.month for d in dates_single]
    d_day_single = [d.day for d in dates_single]

    for plot_date in plot_dates:

        #pdb.set_trace()

        date_match_idx = np.where((np.array(d_year)==plot_date.year) & (np.array(d_month)==plot_date.month) & (np.array(d_day)==plot_date.day))[0]
        #date_match_idx_single = np.where((np.array(d_year_single)==plot_date.year) & (np.array(d_month_single)==plot_date.month) & (np.array(d_day_single)==plot_date.day))[0]
        

        for ds, d in enumerate (date_match_idx):
            try:           
                #plot_dates = dates[d]

                #ds_idx = date_match_idx_single[ds] 

                #pdb.set_trace()

                plot_data = data[:, d]

                p=plot_data[0, :]/100
                T=plot_data[1, :]-273.15
                Td=T-plot_data[2, :]
                h=plot_data[15, :]
                da=dates[d]
                #print T
                #print p
                #print Td

                #pdb.set_trace()

                u_wind,v_wind = u_v_winds(plot_data[3, :], plot_data[4, :])
                
                p_wind = p[~np.isnan(u_wind)]

                u_wind = u_wind[~np.isnan(u_wind)]
                v_wind = v_wind[~np.isnan(v_wind)]     
                
        
                # Create a new figure. The dimensions here give a good aspect ratio
                fig = plt.figure(figsize=(10, 8), frameon=False)
                #fig.patch.set_visible(False)
            
  
        
                tephigram_plot_height=0.85
                tephigram_plot_bottom=.085
                ax = fig.add_axes([.085,tephigram_plot_bottom,.65,tephigram_plot_height], projection='skewx', frameon=False, axisbg='w')
                ax.set_yscale('log')
                plt.grid(True)

                wbax = fig.add_axes([0.75,tephigram_plot_bottom,0.12,tephigram_plot_height],frameon=False, sharey=ax, label='barbs')
                ax_text_box = fig.add_axes([0.85,0.085,.12,tephigram_plot_height], frameon=False, axisbg='w')

                #pdb.set_trace()
                #tmax=math.ceil(nanmax(T)/10)*10
                #tmin=math.floor(nanmin(Td[p>400])/10)*10
                pmax=math.ceil(nanmax(p)/50)*50

                tmax=40.
                tmin=-40.
                
                P=linspace(pmax,pmin,37)

                w = array([0.0001,0.0004,0.001, 0.002, 0.004, 0.007, 0.01, 0.016, 0.024, 0.032, 0.064, 0.128])
                ax.add_mixratio_isopleths(w,P,color='m',ls='-',alpha=.5,lw=0.5)
                ax.add_dry_adiabats(linspace(250,440,20)-273.15,P,color='g',ls='-',alpha=.5,lw=0.8)
                ax.add_moist_adiabats(linspace(8,32,7),P,color='b',ls='-',alpha=.5,lw=0.8)
                ax.other_housekeeping(pmax, pmin, tmax,tmin) 

                # Plot the data using normal plotting functions, in this case using semilogy
 
                ax.semilogy(T[~np.isnan(T)], p[~np.isnan(T)], 'k', linewidth=2)
                ax.semilogy(Td[~np.isnan(Td)], p[~np.isnan(Td)], 'k',linewidth=2)

                row_labels=(
                        'SLAT',
                        'SLON',
                        'SELV',
                        'SHOW',
                        'LIFT',
                        'LFTV',
                        'SWET',
                        'KINX',
                        'CTOT',
                        'VTOT',
                        'TOTL',
                        'CAPE',
                        'CINS',
                        'CAPV',
                        'CINV',
                        'LFCT',
                        'LFCV',
                        'BRCH',
                        'BRCV',
                        'LCLT',
                        'LCLP',
                        'MLTH', 
                        'MLMR',
                        'THCK',
                        'PWAT')

                variable='pbl_pressure'
                var_index = variable_name_index_match(variable, variable_list_line)
                print load_file['single_value_vars'].shape 

                
                #pbl_pressure =  load_file['single_value_vars'][var_index, ds_idx]

                #print pbl_pressure

                # EQLV, pp, lclp,lfcp, lclt, delta_z, CAPE, CIN=CapeCinPBLInput(p, T, Td, h, st_height, pbl_pressure/100)

                # print lclp


                # table_vals=(
                # #'%s' % station_name,
                # #'Climatology - Week beg. %s' % da, 
                # '%s' % la,
                # '%s' % lo,
                # '%s' % st_height,
                # '%.1f' % ShowalterIndex(T, Td, p),                       #      ['Showalter index',
                # '%.1f' % LiftedIndex(T, Td, p, h, st_height),          #     'Lifted index', 
                # '--',                                                                    #    'LIFT computed using virtual temperature',
                # '--',                                                                    #   'SWEAT index',
                # '%.1f' % KIndex(T, Td, p),          #    'K index',
                # '%.1f' % CrossTotalsIndex(T, Td, p),                   #   'Cross totals index',
                # '%.1f' % VerticalTotalsIndex(T, p),                      #   'Vertical totals index',
                # '%.1f' % TotalTotalsIndex(T, Td, p),                     #   'Total totals index',
                # '%.1f' % CAPE,                                                   #   'CAPE',
                # '%.1f' % CIN,                                                      #   'CIN',
                # '--',                                                                    #   'CAPE using virtual temperature',
                # '--',                                                                    #   'CINS using virtual temperature',
                # '%.1f' % lfcp,                                                     #   'Level of free convection',
                # '--',                                                                    #   'LFCT using virtual temperature',
                # '--' ,                                                                   #   'Bulk Richardson number',
                # '--',                                                                    #   'Bulk richardson using CAPV',
                # '%.1f' % lclt,                                                        #  'Temp [K] of the Lifted Condensation Level',
                # '%.1f' % lclp ,                                                       #  'Pres [hPa] of the Lifted Condensation Level',
                # '--',                                                                    #  'Mean mixed layer potential temperature', 
                # '--',                                                                    #  'Mean mixed layer mixing ratio',
                # '--',                                                                    #  '1000 hPa to 500 hPa thickness',
                # '--')                                                                    #  'Precipitable water [mm] for entire sounding'] 
            
                # Wind barbs

                #pdb.set_trace()

                barbs_idx=np.logspace(np.log10(10),np.log10(max(len(u_wind))),num=32).astype(int)
                #barbs_idx=np.logspace(np.log10(10),np.log10(max(len(u_wind))),num=32).astype(int)
                wbax.set_yscale('log')
                wbax.xaxis.set_ticks([],[])
                wbax.yaxis.grid(True,ls='-',color='y',lw=0.5)

                wbax.set_xlim(-1.5,1.5)
                wbax.get_yaxis().set_visible(False)
                wbax.set_ylim(pmax+100,pmin)
                
                #pdb.set_trace()

                wbax.barbs((zeros(p_wind.shape)),p_wind, u_wind, v_wind)

                # Disables the log-formatting that comes with semilogy
                ax.yaxis.set_major_formatter(ScalarFormatter())
                ax.set_yticks(linspace(100,1000,10))
                ax.set_ylim(pmax,pmin)

                ax.set_xlim(tmin,tmax)
                ax.xaxis.set_ticks([],[])
     
                ax_text_box.xaxis.set_visible(False)
                ax_text_box.yaxis.set_visible(False)
                for tick in wbax.yaxis.get_major_ticks():
                    # tick.label1On = False
                    pass
                    #wbax.get_yaxis().set_tick_params(size=0,color='y')
                

                y_loc=1.

                max_string_length = max([len(line) for line in row_labels])
        
                # for t,r in zip(row_labels,table_vals):
                #     label_rightjust=('{:>%i}' % max_string_length).format(t)
      
                #     ax_text_box.text(0.5, y_loc, ' %s:' % (label_rightjust), size=8, horizontalalignment='right')
                #     ax_text_box.text(0.5, y_loc, ' %s' % (r), size=8, horizontalalignment='left')
                #     y_loc-=0.04

                fig.text(.02,0.965, '%s  %s' %(stat, station_name), size=12, horizontalalignment='left')
                fig.text(.02,0.035, '%s ' %(da.strftime('%Y-%m-%d %H:%M')), size=12, horizontalalignment='left')
                #plt.show()
                plt.savefig('/nfs/a90/eepdw/Figures/Radiosonde/Tephigrams/%s_%s_%s_Skew_T_Vars_To_Right_Barbs.png' % (station_name.replace('/','_').replace(' ', '_'), stat, dates[d].strftime('%Y%m%d')))
                plt.close()
            except Exception:
                print PrintException()
