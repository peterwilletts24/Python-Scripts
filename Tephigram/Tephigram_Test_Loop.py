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
imp.load_source('TephigramPlot', '/nfs/see-fs-01_users/eepdw/python_scripts/Tephigram/Tephigram_Test.py')
from TephigramPlot import *
from SoundingRoutines import *

pmax=1000.
pmin=200.

station_list_cs=[43003, 43014, 42867, 43371, 43353, 43285,  43192, 43150, 42339, 40990, 40948]

date_min=datetime.datetime(1960,5,1,0,0,0)
date_max=datetime.datetime(2014,10,1,0,0,0)

delta = relativedelta(weeks=+1)    
        
variable_list={'pressures': 0, 'temps':1, 'dewpoints':2, 'winddirs':3, 'windspeeds':4, 'pot_temp':5, 
               'sat_vap_pres':6, 'vap_press':7, 'rel_hum':8, 'wvmr':9, 'sp_hum':10, 'sat_temp':11, 
               'theta_e':12, 'theta_e_sat':13, 'theta_e_minus_theta_e_sat':14}

station_list_search='/nfs/a90/eepdw/Data/Observations/Radiosonde_downloaded_from_NOAA_GUAN/igra-stations.txt'

station_metadata=[]
f = open(station_list_search,'r')
for line in f:
     line = line.strip()
     line=re.sub(r'([A-Z])\s([A-Z])', r'\1_\2',line)
     line=re.sub(r'([A-Z])\s\s([A-Z])', r'\1_\2',line)
     station_metadata.append(line.split())
f.close()

    # Parse the data
for stat in station_list_cs:

    station_name,la,lo, st_height=station_info_search(stat)
    
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

            u_wind,v_wind = u_v_winds(data[bin,3,:], data[bin,4,:])
        
            # Create a new figure. The dimensions here give a good aspect ratio
            fig = plt.figure(figsize=(7, 16), frameon=False)
            #fig.patch.set_visible(False)
            
            tephigram_plot_height=0.47
            tephigram_plot_bottom=0.5
            ax = fig.add_axes([.085,tephigram_plot_bottom,.75,tephigram_plot_height], projection='skewx', frameon=False, axisbg='w')
            ax.set_yscale('log')
            plt.grid(True)

            #pdb.set_trace()
            tmax=math.ceil(nanmax(T)/10)*10
            tmin=math.floor(nanmin(Td[p>400])/10)*10
    
            P=linspace(pmax,pmin,37)

            w = array([0.0001,0.0004,0.001, 0.002, 0.004, 0.007, 0.01, 0.016, 0.024, 0.032, 0.064, 0.128])
            ax.add_mixratio_isopleths(w,P,color='m',ls='-',alpha=.5,lw=0.5)
            ax.add_dry_adiabats(linspace(250,440,20)-273.15,P,color='g',ls='-',alpha=.5,lw=0.8)
            ax.add_moist_adiabats(linspace(8,32,7),P,color='b',ls='-',alpha=.5,lw=0.8)
            ax.other_housekeeping(pmax, pmin, tmax,tmin) 

            # Plot the data using normal plotting functions, in this case using semilogy
 
            ax.semilogy(T, p, 'k', linewidth=2)
            ax.semilogy(Td, p, 'k',linewidth=2)

            row_labels=('Station name',
                        'Date:',
                        'Latitude',
                        'Longitude',
                        'Elevation(m)',
                        'Showalter index',
                        'Lifted index',
                        'LIFT computed using virtual temperature',
                        'SWEAT index',
                        'K index',
                        'Cross totals index',
                        'Vertical totals index',
                        'Total totals index',
                        'CAPE',
                        'CIN',
                        'CAPE using virtual temperature',
                        'CINS using virtual temperature',
                        'Level of free convection',
                        'LFCT using virtual temperature',
                        'Bulk Richardson number',
                        'Bulk richardson using CAPV',
                        'Temp [K] of the Lifted Condensation Level',
                        'Pres [hPa] of the Lifted Condensation Level',
                        'Mean mixed layer potential temperature', 
                        'Mean mixed layer mixing ratio',
                        '1000 hPa to 500 hPa thickness',
                        'Precipitable water [mm] for entire sounding')

            EQLV, pp, lclp,lfcp, lclt, delta_z, CAPE, CIN=CapeCin(p, T, Td, h, st_height)

            table_vals=(
                '%s' % station_name,
                'Climatology - Week beg. %s' % da, 
                '%s' % la,
                '%s' % lo,
                '%s' % st_height,
                '%s' % ShowalterIndex(T, Td, p),                       #      ['Showalter index',
                '%s' % LiftedIndex(T, Td, p, h, st_height),          #     'Lifted index', 
                '--',                                                                    #    'LIFT computed using virtual temperature',
                '--',                                                                    #   'SWEAT index',
                '%s' % LiftedIndex(T, Td, p, h, st_height),          #    'K index',
                '%s' % CrossTotalsIndex(T, Td, p),                   #   'Cross totals index',
                '%s' % VerticalTotalsIndex(T, p),                      #   'Vertical totals index',
                '%s' % TotalTotalsIndex(T, Td, p),                     #   'Total totals index',
                '%s' % CAPE,                                                   #   'CAPE',
                '%s' % CIN,                                                      #   'CIN',
                '--',                                                                    #   'CAPE using virtual temperature',
                '--',                                                                    #   'CINS using virtual temperature',
                '%s' % lfcp,                                                     #   'Level of free convection',
                '--',                                                                    #   'LFCT using virtual temperature',
                '--' ,                                                                   #   'Bulk Richardson number',
                '--',                                                                    #   'Bulk richardson using CAPV',
                '%s' % lclt,                                                        #  'Temp [K] of the Lifted Condensation Level',
                '%s' % lclp ,                                                       #  'Pres [hPa] of the Lifted Condensation Level',
                '--',                                                                    #  'Mean mixed layer potential temperature', 
                '--',                                                                    #  'Mean mixed layer mixing ratio',
                '--',                                                                    #  '1000 hPa to 500 hPa thickness',
                '--')                                                                    #  'Precipitable water [mm] for entire sounding'] 


            skip=max(len(u_wind)/32)
            bloc=0.5

            # Wind barbs

            wbax=fig.add_axes([0.88,tephigram_plot_bottom,0.12,tephigram_plot_height],frameon=False, sharey=ax, label='barbs')

            wbax.xaxis.set_ticks([],[])
            wbax.yaxis.grid(True,ls='-',color='y',lw=0.5)
            for tick in wbax.yaxis.get_major_ticks():
                # tick.label1On = False
                pass
                #wbax.get_yaxis().set_tick_params(size=0,color='y')
                wbax.set_xlim(-1.5,1.5)
                wbax.get_yaxis().set_visible(False)
                wbax.set_ylim(pmax,pmin)
                wbax.set_yscale('log')
                wbax.barbs((zeros(p.shape))[::skip],p[::skip], u_wind[::skip], v_wind[::skip])

                # Disables the log-formatting that comes with semilogy
                ax.yaxis.set_major_formatter(ScalarFormatter())
                ax.set_yticks(linspace(100,1000,10))
                ax.set_ylim(pmax,pmin)

                ax.set_xlim(tmin,tmax)
                ax.xaxis.set_ticks([],[])

                ax_text_box = fig.add_axes([.085,0.1,.75,.35], frameon=False, axisbg='w')
     
                ax_text_box.xaxis.set_visible(False)
                ax_text_box.yaxis.set_visible(False)

                y_loc=1.

                max_string_length = max([len(line) for line in row_labels])
        
                for t,r in zip(row_labels,table_vals):
                    label_rightjust=('{:>%i}' % max_string_length).format(t)
      
                    ax_text_box.text(0.5, y_loc, ' %s:' % (label_rightjust), size=8, horizontalalignment='right')
                    ax_text_box.text(0.5, y_loc, ' %s' % (r), size=8, horizontalalignment='left')
                    y_loc-=0.04

                #plt.show()
                plt.savefig('/nfs/a90/eepdw/Figures/Radiosonde/Tephigrams/Weekly_Climatology/%s_%s_%s_Skew_T.png' % (station_name.replace('/','_').replace(' ', '_'), stat, da.strftime('%Y%m%d')))
        except Exception:
            print PrintException()
