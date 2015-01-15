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

imp.load_source('IrisFunctions', '/nfs/see-fs-01_users/eepdw/python_scripts/modules/IrisFunctions.py')
from IrisFunctions import *

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

top_dir='/nfs/a90/eepdw/Data/EMBRACE'
save_dir='/nfs/a90/eepdw/Figures/EMBRACE/Diurnal'

pp_file = '_408_on_p_levs_mean_by_hour_diurnal_'

region = 'monsoon_trough'

types_of_plot=['large_domain_only', '8_and_12_km_only', 'all']

#types_of_plot=['8_and_12_km_only', 'all']
types_of_plot=['all']

# lon_max = 101.866 
# lon_min = 64.115
# lat_max= 33.
# lat_min=-6.79

#trmm_dir = '/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/Diurnal/'
#trmm_file = "trmm_diurnal_average_lat_%s_%s_lon_%s_%s.npz" % (lat_min,lat_max, lon_min, lon_max)

#############



# Make own time x-axis (UTC)

#d = matplotlib.dates.drange(datetime.datetime(2011, 8, 21, 6,30), datetime.datetime(2011, 8, 22, 6, 30), timedelta(hours=1))

formatter = matplotlib.dates.DateFormatter('%H:%M')

def main():
  for type_of_plot in types_of_plot:
    if type_of_plot=='large_domain_only':
        experiment_ids_p = [ 'djznw', 'dkjxq', 'djznq' ] # Params
        experiment_ids_e = ['dkhgu', 'dkbhu'] # Explicit
    if type_of_plot=='8_and_12_km_only':
        experiment_ids_p = [ 'djznw' 'dkmbq', 'dklzq' ] # Params
        experiment_ids_e = ['dklyu', 'dklwu', 'dkbhu'] # Explicit
    if type_of_plot=='all':
        experiment_ids_p = ['djznw', 'djzny', 'djznq', 'dklzq', 'dkmbq', 'dkjxq' ] # Most of Params
        experiment_ids_e = ['dklwu', 'dklyu', 'djzns', 'dkbhu', 'djznu', 'dkhgu'] # Most of Explicit

    NUM_COLOURS = 16
    cmap=cm.get_cmap(cm.Set1, NUM_COLOURS)
 #cgen = (cmap(1.*i/NUM_COLORS) for i in range(NUM_COLORS))

    for ls in ['land']:
 # for ls in ['sea']:
     fig = plt.figure(figsize=(12,6))
     ax = fig.add_subplot(111)
     #legendEntries=[]
     #legendtext=[]

     #if ls=='land':
      #bbox_anchor=-0.25
     #else:
     
     bbox_anchor=0.1


       # Change the legend label colors to almost black
     #texts = l0.texts
     #for t in texts:
      #t.set_color('#262626')

     legendEntries=[]
     legendtext=[] 


     for c, experiment_id in enumerate(experiment_ids_p):

         expmin1 = experiment_id[:-1]
  

         if (experiment_id=='djznw'):
          print experiment_id
          #colour = cmap(1.*1/NUM_COLOURS)
          colour = '#262626'
          linewidth=0.5
          linestylez='--'
         if (experiment_id=='djzny'):
          print experiment_id
          colour = cmap(1.*3/NUM_COLOURS)
          linewidth=0.5
          linestylez='--'
         if ((experiment_id=='djznq') or (experiment_id=='dkjxq')):
          print experiment_id
          colour = cmap(1.*11/NUM_COLOURS)
          linewidth=0.8
          if (experiment_id=='djznq'):
              linestylez='--'
          if (experiment_id=='dkjxq'):
              linestylez=':'
              
         if ((experiment_id=='dklzq') or (experiment_id=='dklwu')):
          print experiment_id
          colour = cmap(1.*7/NUM_COLOURS)
          linewidth=1
          if (experiment_id=='dklzq'):
              linestylez='--'
          if (experiment_id=='dklwu'):
              linestylez='-'
         if ((experiment_id=='dklyu') or (experiment_id=='dkmbq')):
          print experiment_id
          colour = cmap(1.*9/NUM_COLOURS)
          linewidth=1.3
          if (experiment_id=='dkmbq'):
              linestylez='--'
          if (experiment_id=='dklyu'):
              linestylez='-'
         if (experiment_id=='djzns'):
          print experiment_id
          colour = cmap(1.*11/NUM_COLOURS)
          linewidth=1.6
          linestylez='-'
         if ((experiment_id=='dkbhu')or (experiment_id=='dkhgu')):
          print experiment_id
          colour = cmap(1.*13/NUM_COLOURS)
          linewidth=1.9
          if (experiment_id=='dkbhu'):
              linestylez='-'
          if (experiment_id=='dkhgu'):
              linestylez=':'
         if (experiment_id=='djznu'):
          print experiment_id
          colour = cmap(1.*15/NUM_COLOURS)
          linewidth=2.
          linestylez='-'
         try:
          plotnp = np.load('%s/%s/%s/%s%s%s.npz' % (top_dir, expmin1, experiment_id, experiment_id, pp_file, region))

          #if (ls != 'total'):
                # Make own time x-axis (local)

          #pdb.set_trace()

          hours_to_plot = ConvertHoursSince1970ToDatetime(plotnp['time'])+utc_to_local
          hours_to_plot = np.array([datetime.datetime(2011, 8, 19, h.hour, h.minute) for h in hours_to_plot])

          hour_arg_sort=np.argsort([h.hour for h in hours_to_plot])

             #time_sort = plotnp[1][hour_arg_sort]
          data_sort =  plotnp['data'][:,-3][hour_arg_sort]

          l, = plt.plot_date(hours_to_plot[hour_arg_sort], data_sort, label=model_name_convert_legend.main(experiment_id), linewidth=linewidth, linestyle=linestylez, marker='', markersize=2, fmt='', color=colour)
          #else:
          # l, = plt.plot_date(d, plotnp*3600, label=model_name_convert_legend.main(experiment_id), linewidth=linewidth, linestyle=linestylez, marker='', markersize=2, fmt='', color=colour) 

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


            if (experiment_id=='djznw'):
             print experiment_id
             #colour = cmap(1.*1/NUM_COLOURS)
             colour = '#262626'
             linewidth=0.5
             linestylez='--'
            if (experiment_id=='djzny'):
             print experiment_id
             colour = cmap(1.*3/NUM_COLOURS)
             linewidth=0.5
             linestylez='--'
            if ((experiment_id=='djznq') or (experiment_id=='dkjxq')):
             print experiment_id
             colour = cmap(1.*11/NUM_COLOURS)
             linewidth=0.8
             if (experiment_id=='djznq'):
              linestylez='--'
             if (experiment_id=='dkjxq'):
              linestylez=':'
              
            if ((experiment_id=='dklzq') or (experiment_id=='dklwu')):
             print experiment_id
             colour = cmap(1.*7/NUM_COLOURS)
             linewidth=1
            if (experiment_id=='dklzq'):
             linestylez='--'
            if (experiment_id=='dklwu'):
             linestylez='-'
            if ((experiment_id=='dklyu') or (experiment_id=='dkmbq')):
             print experiment_id
             colour = cmap(1.*9/NUM_COLOURS)
             linewidth=1.3
            if (experiment_id=='dkmbq'):
             linestylez='--'
            if (experiment_id=='dklyu'):
             linestylez='-'
            if (experiment_id=='djzns'):
             print experiment_id
             colour = cmap(1.*11/NUM_COLOURS)
             linewidth=1.6
             linestylez='-'
            if ((experiment_id=='dkbhu')or (experiment_id=='dkhgu')):
             print experiment_id
             colour = cmap(1.*13/NUM_COLOURS)
             linewidth=1.9
             if (experiment_id=='dkbhu'):
              linestylez='-'
             if (experiment_id=='dkhgu'):
              linestylez=':'
            if (experiment_id=='djznu'):
             print experiment_id
             colour = cmap(1.*15/NUM_COLOURS)
             linewidth=2.
             linestylez='-'


            expmin1 = experiment_id[:-1]

            try:
              plotnp = np.load('%s/%s/%s/%s%s%s.npz' % (top_dir, expmin1, experiment_id, experiment_id, pp_file, region))

              #if (ls != 'total'):
                # Make own time x-axis (local)

              #pdb.set_trace()

              hours_to_plot = ConvertHoursSince1970ToDatetime(plotnp['time'])+utc_to_local
              hours_to_plot = np.array([datetime.datetime(2011, 8, 19, h.hour, h.minute) for h in hours_to_plot])

              hour_arg_sort=np.argsort([h.hour for h in hours_to_plot])

              #time_sort = plotnp[1][hour_arg_sort]
              data_sort =  plotnp['data'][:,-3][hour_arg_sort]

              l, = plt.plot_date(hours_to_plot[hour_arg_sort], data_sort, label=model_name_convert_legend.main(experiment_id), linewidth=linewidth, linestyle=linestylez, marker='', markersize=2, fmt='', color=colour)
              #else:try:
             
              plotnp = np.load('%s/%s/%s/%s%s%s.npz' % (top_dir, expmin1, experiment_id, experiment_id, pp_file, region))

      #plotnp = np.ort(pnp, axis=1)
             #if (ls != 'total'):

             # Make own time x-axis (local)

          
              legendEntries.append(l)
              legendtext.append('%s' % (model_name_convert_legend.main(experiment_id)))

            except Exception, e:
             print e
             pass

     l2=plt.legend(legendEntries, legendtext, title='Explicit', loc=9, frameon=False, prop={'size':12}, 
                   bbox_to_anchor=(0.155+bbox_anchor, 0,1, 1))
     plt.gca().add_artist(l1)
     #plt.gca().add_artist(l0)
     plt.gca().xaxis.set_major_formatter(formatter)

             # Change the legend label colors to almost black
     texts = l2.texts
     for t in texts:
      t.set_color('#262626')

     plt.xlabel('Time (Local)')
     plt.ylabel('${Height (m)}$')

     if pp_file == '1201':
       diag_title = 'Net Down Surface SW Flux'
     if pp_file == '2201':
       diag_title = 'Net Down Surface LW Flux'
     if pp_file == '3217':
       diag_title = 'Surface Heat Flux'
     if pp_file == '3234':
       diag_title = 'Surface Latent Heat Flux'
     if pp_file == '_408_on_p_levs_mean_by_hour_diurnal_':
       diag_title = '925 hPa Height'


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

     plt.savefig('%s/EMBRACE_Diurnal_%s_%s_%s_%s_notitle.png' % (save_dir, pp_filenodot, ls, type_of_plot, region),
                 format='png', bbox_inches='tight')
                  


     plt.title('\n'.join(wrap('%s' % (t.title()), 1000,replace_whitespace=False)), fontsize=16)
      #plt.show()
  
     plt.savefig('%s/EMBRACE_Diurnal_%s_%s_%s_%s.png' % (save_dir, pp_filenodot, ls, type_of_plot, region),
                 format='png', bbox_inches='tight')
     plt.close()

     #except Exception, e:
      #print e


if __name__ == '__main__':
   main()
       



