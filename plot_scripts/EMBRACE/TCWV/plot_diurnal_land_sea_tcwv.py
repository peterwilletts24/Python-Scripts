"""

Load npy xy, plot and save


"""

import os, sys

import matplotlib

 #matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!

import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm

from matplotlib import rc
from matplotlib.font_manager import FontProperties
from matplotlib import rcParams
from matplotlib import cm

import pdb

import imp
imp.load_source('PlotFunctions', '/nfs/see-fs-01_users/eepdw/python_scripts/modules/PlotFunctions.py')
#from IrisFunctions import *
from PlotFunctions import *

rc('text', usetex=True)

rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True

rc('font', family = 'serif', serif = 'cmr10')

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 14}

matplotlib.rc('font', **font)

import numpy as np

from datetime import timedelta
import datetime

import math

import imp

import re
from textwrap import wrap

import iris.analysis.geometry
from shapely.geometry import Polygon

model_name_convert_legend = imp.load_source('util', '/nfs/see-fs-01_users/eepdw/python_scripts/modules/model_name_convert_legend.py')
#unrotate = imp.load_source('util', '/home/pwille/python_scripts/modules/unrotate_pole.py')
utc_to_local=datetime.timedelta(hours=5, minutes=30)
###############
# Things to change

time_difference = 5.5

pp_file = 'TCWV'

top_dir='/nfs/a90/eepdw/Data/EMBRACE'
save_dir='/nfs/a90/eepdw/Figures/EMBRACE/Diurnal/%s' % pp_file

types_of_plot=['large_domain_only', '8_and_12_km_only', 'all']

#types_of_plot=['8_and_12_km_only', 'all']
#types_of_plot=['all']

# lon_max = 101.866 
# lon_min = 64.115
# lat_max= 33.
# lat_min=-6.79

#trmm_dir = '/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/Diurnal/'
#trmm_file = "trmm_diurnal_average_lat_%s_%s_lon_%s_%s.npz" % (lat_min,lat_max, lon_min, lon_max)

#############



# Make own time x-axis (UTC)

#d = matplotlib.dates.drange(datetime.datetime(2011, 8, 21, 6,30), datetime.datetime(2011, 8, 22, 6, 30), timedelta(hours=1))

#formatter = matplotlib.dates.DateFormatter('%H:%M')

def main():
  for type_of_plot in types_of_plot:
    if type_of_plot=='large_domain_only':
        experiment_ids_p = [ 'djznw', 'dkjxq', 'djznq' ] # Params
        experiment_ids_e = ['dkhgu', 'dkbhu'] # Explicit
    if type_of_plot=='8_and_12_km_only':
        experiment_ids_p = [ 'djznw', 'dkmbq', 'dklzq' ] # Params
        experiment_ids_e = ['dklyu', 'dklwu', 'dkbhu'] # Explicit
    if type_of_plot=='all':
        experiment_ids_p = ['djznw', 'djzny', 'djznq', 'dklzq', 'dkmbq', 'dkjxq' ] # Most of Params
        experiment_ids_e = ['dklwu', 'dklyu', 'djzns', 'dkbhu', 'djznu', 'dkhgu'] # Most of Explicit

    for ls in ['land', 'sea', 'total']:
 # for ls in ['sea']:
     fig = plt.figure(figsize=(12,6))
     ax = fig.add_subplot(111)
     #legendEntries=[]
     #legendtext=[]

     if ls=='sea':
      bbox_anchor=-0.25
     else:
      bbox_anchor=0


       # Change the legend label colors to almost black
     #texts = l0.texts
     #for t in texts:
      #t.set_color('#262626')

     legendEntries=[]
     legendtext=[] 


     for c, experiment_id in enumerate(experiment_ids_p):

         expmin1 = experiment_id[:-1]
  
         try:

                colour, linewidth, linestylez = LinePlotEMBRACEExperimentID(experiment_id)

         except Exception:
                 print ' colour not assigned %s' % experiment_id



         try:
          plot_data = np.load('%s/%s/%s/%s_%s_%s_hourly_mean.npz' % (top_dir, expmin1, experiment_id, pp_file, experiment_id, ls))

         

          hours_local_ascend = np.where(plot_data['hours']+time_difference<24, plot_data['hours']+time_difference, plot_data['hours']+time_difference-24)
          x_sort=np.argsort(hours_local_ascend)

          l, = plt.plot(hours_local_ascend[x_sort], -plot_data['means_hourly'][x_sort]*1000,  linewidth=linewidth, linestyle=linestylez, marker='', markersize=2, color=colour)

          legendEntries.append(l)
          legendtext.append('%s' % (model_name_convert_legend.main(experiment_id)))
  

         except Exception, e:
          print e
          pass

     l1=plt.legend(legendEntries, legendtext, title='Parametrised', loc=9, frameon=False, prop={'size':12}, 
                 bbox_to_anchor=(0+bbox_anchor, 0,1, 1))

     # Change the legend label colors to almost black
     texts = l1.texts
     for t in texts:
      t.set_color('#262626')

      legendEntries=[]
     legendtext=[]
 
     c1=0
     for c, experiment_id in enumerate(experiment_ids_e):

            expmin1 = experiment_id[:-1]

            try:

                colour, linewidth, linestylez = LinePlotEMBRACEExperimentID(experiment_id)

            except Exception:
                 print experiment_id

            try:
   
                    plot_data = np.load('%s/%s/%s/%s_%s_%s_hourly_mean.npz' % (top_dir, expmin1, experiment_id, pp_file, experiment_id, ls))

                    hours_local_ascend = np.where(plot_data['hours']+time_difference<24, plot_data['hours']+time_difference, plot_data['hours']+time_difference-24)
                    x_sort=np.argsort(hours_local_ascend)

                    l,=plt.plot(hours_local_ascend[x_sort], -plot_data['means_hourly'][x_sort]*1000,  linewidth=linewidth, linestyle=linestylez, marker='', markersize=2, color=colour)

                    legendEntries.append(l)
                    legendtext.append('%s' % (model_name_convert_legend.main(experiment_id)))

            except Exception, e:
             print e
             pass

     l2=plt.legend(legendEntries, legendtext, title='Explicit', loc=9, frameon=False, prop={'size':12}, 
                   bbox_to_anchor=(0.155+bbox_anchor, 0,1, 1))
     plt.gca().add_artist(l1)
     #plt.gca().add_artist(l0)
     #plt.gca().xaxis.set_major_formatter(formatter)

             # Change the legend label colors to almost black
     texts = l2.texts
     for t in texts:
      t.set_color('#262626')

     plt.xlabel('Time (Local)')
     #plt.ylabel('${W m^-2}$')
     plt.ylabel('${mm}$') 

     if pp_file == '1201':
       diag_title = 'Net Down Surface SW Flux'
     if pp_file == '2201':
       diag_title = 'Net Down Surface LW Flux'
     if pp_file == '3217':
       diag_title = 'Surface Heat Flux'
     if pp_file == '3234':
       diag_title = 'Surface Latent Heat Flux'
     if pp_file == 'TCWV':
       diag_title = 'Total Column Water Vapour'


     title='Domain Averaged   %s - %s' % (diag_title, ls)

     t=re.sub('(.{68} )', '\\1\n', str(title), 0, re.DOTALL)
     t = re.sub(r'[(\']', ' ', t)
     t = re.sub(r'[\',)]', ' ', t)
            
     pp_filenodot= pp_file.replace(".", "")

     # Bit of formatting 
     # Set colour of axis lines
     spines_to_keep = ['bottom', 'left']
     for spine in spines_to_keep:
      ax.spines[spine].set_linewidth(0.5)
      ax.spines[spine].set_color('#262626')
      # Remove top and right axes lines ("spines")
      spines_to_remove = ['top', 'right']
     for spine in spines_to_remove:
       ax.spines[spine].set_visible(False)
       # Get rid of ticks. The position of the numbers is informative enough of
       # the position of the value.
       ax.xaxis.set_ticks_position('none')
       ax.yaxis.set_ticks_position('none')
             
      # Change the labels to the off-black
     ax.xaxis.label.set_color('#262626')
     ax.yaxis.label.set_color('#262626')

     if not os.path.exists(save_dir): os.makedirs(save_dir)

     #plt.savefig('%s/EMBRACE_Diurnal_%s_%s_%s_notitle.png' % (save_dir, pp_filenodot, ls, type_of_plot),
         #        format='png', bbox_inches='tight')
                  
     plt.xlim(0,24)

    
     plt.title('%s %s Diurnal Mean' % (pp_filenodot, ls))

     plt.savefig('%s/EMBRACE_Diurnal_%s_%s_%s.png' % (save_dir, pp_filenodot, ls, type_of_plot),
                format='png', bbox_inches='tight')
     plt.close()

     #except Exception, e:
      #print e


if __name__ == '__main__':
   main()
       



