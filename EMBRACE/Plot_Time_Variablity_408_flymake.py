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

import linecache

import datetime

def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print 'EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj)

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

imp.load_source('IrisFunctions', '/nfs/see-fs-01_users/eepdw/python_scripts/modules/IrisFunctions.py')
from IrisFunctions import *

#unrotate = imp.load_source('util', '/home/pwille/python_scripts/modules/unrotate_pole.py')
utc_to_local=datetime.timedelta(hours=5, minutes=30)
###############
# Things to change

top_dir='/nfs/a90/eepdw/Data/EMBRACE'
save_dir='/nfs/a90/eepdw/Figures/EMBRACE/TimeVariability'

file_name = 'TimeVar_408_land_domain_constrain__and_oro_not_greater_than_data_monsoon_trough.npz'
    
diag_title = 'Area Averaged Rainfall'
area='monsoon_trough'

types_of_plot=['large_domain_only', '8_and_12_km_only', 'all']

#types_of_plot=['8_and_12_km_only', 'all']
types_of_plot=['all']
types_of_plot=['large_domain_only', '8_and_12_km_only', 'all', '8_and_12_km_plus', '8_and_12_km_par_only', '8_and_12_km_exp_only']
formatter = matplotlib.dates.DateFormatter('%d %b')

# lon_max = 101.866 
# lon_min = 64.115
# lat_max= 33.
# lat_min=-6.79

#trmm_dir = '/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/Diurnal/'
#trmm_file = "trmm_diurnal_average_lat_%s_%s_lon_%s_%s.npz" % (lat_min,lat_max, lon_min, lon_max)


def main():
  for type_of_plot in types_of_plot:
    if type_of_plot=='large_domain_only':
        experiment_ids_p = [ 'djznw', 'dkjxq', 'djznq' ] # Params
        experiment_ids_e = ['dkhgu', 'dkbhu'] # Explicit

    if type_of_plot=='8_and_12_km_par_only':
        experiment_ids_p = [ 'dkmbq', 'dklzq' ] # Params
        experiment_ids_e = [] # Explicit        

    if type_of_plot=='8_and_12_km_exp_only':
        experiment_ids_p = [] # Params
        experiment_ids_e = ['dklyu', 'dklwu', 'dkmgw'] # Explicit  
    if type_of_plot=='8_and_12_km_plus':
        experiment_ids_p = [ 'djznw' ,'dkmbq', 'dklzq' ] # Params
        experiment_ids_e = ['dklyu', 'dkmgw', 'dklwu', 'dkbhu'] # Explicit

    if type_of_plot=='8_and_12_km_only':
        experiment_ids_p = [ 'djznw', 'dkmbq', 'dklzq' ] # Params
        experiment_ids_e = ['dklyu', 'dkmgw', 'dklwu', 'dkbhu'] # Explicit
    if type_of_plot=='all':
        experiment_ids_p = ['djznw', 'djzny', 'djznq', 'dklzq', 'dkmbq', 'dkjxq' ] # Most of Params
        experiment_ids_e = ['dklwu', 'dklyu', 'dkmgw', 'djzns', 'dkbhu', 'djznu', 'dkhgu'] # Most of Explicit

    NUM_COLOURS = 16
    cmap=cm.get_cmap(cm.Set1, NUM_COLOURS)
 #cgen = (cmap(1.*i/NUM_COLORS) for i in range(NUM_COLORS))

    #for ls in ['land', 'sea', 'total']:
    for ls in ['land']:
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
  

         if (experiment_id=='djznw'):
          print experiment_id
          colour = cmap(1.*1/NUM_COLOURS)
          linewidth=0.2
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
          plotnp = np.load('%s/%s/%s/%s' % (top_dir, expmin1, experiment_id, file_name))

          #if (ls != 'total'):
                # Make own time x-axis (local)

          #pdb.set_trace()

          time_arg_sort=np.argsort(plotnp['time_coords'])

             #time_sort = plotnp[1][hour_arg_sort]
          try:
                 data_sort =  plotnp['data'][:,-3][time_arg_sort]
          except Exception:
                 data_sort =  plotnp['data'][time_arg_sort]

          #minute_min,hour_min = math.modf(plotnp['time_coords'].points.min()) 
          #minute_max,hour_max = math.modf(plotnp['time_coords'].points.max()) 

          plot_dates = ConvertHoursSince1970ToDatetime(plotnp['time_coords'][time_arg_sort])

          #pdb.set_trace()

          l, = plt.plot_date(plot_dates, data_sort, label=model_name_convert_legend.main(experiment_id), linewidth=linewidth, linestyle=linestylez, marker='', markersize=2, fmt='', color=colour)
          #else:
          # l, = plt.plot_date(d, plotnp*3600, label=model_name_convert_legend.main(experiment_id), linewidth=linewidth, linestyle=linestylez, marker='', markersize=2, fmt='', color=colour) 

          legendEntries.append(l)
          legendtext.append('%s' % (model_name_convert_legend.main(experiment_id)))
  

         except Exception, e:
          print e
          PrintException()
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
             colour = cmap(1.*1/NUM_COLOURS)
             linewidth=0.2
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
             
             plotnp = np.load('%s/%s/%s/%s' % (top_dir, expmin1, experiment_id, file_name))

             time_arg_sort=np.argsort(plotnp['time_coords'])

             #time_sort = plotnp[1][hour_arg_sort]
             try:
                 data_sort =  plotnp['data'][:,-3][time_arg_sort]
             except Exception:
                 data_sort =  plotnp['data'][time_arg_sort]
             #minute_min,hour_min = math.modf(plotnp['time_coords'].points.min()) 
             #minute_max,hour_max = math.modf(plotnp['time_coords'].points.max()) 

             #pdb.set_trace()

             plot_dates = ConvertHoursSince1970ToDatetime(plotnp['time_coords'][time_arg_sort])

             l, = plt.plot_date(plot_dates, data_sort, label=model_name_convert_legend.main(experiment_id), linewidth=linewidth, linestyle=linestylez, marker='', markersize=2, fmt='', color=colour)
             #else:

      
             
          
             #else:
              #l, = plt.plot_date(d, plotnp*3600, label='%s' % (model_name_convert_legend.main(experiment_id)), linewidth=linewidth, linestyle=linestylez, marker='', markersize=2, fmt='', color=colour)
             legendEntries.append(l)
             legendtext.append('%s' % (model_name_convert_legend.main(experiment_id)))

            except Exception, e:
             print e
             PrintException()
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
     plt.ylabel('${mm h^-1}$')

     title='Domain Averaged   %s - %s' % (diag_title, ls)

     t=re.sub('(.{68} )', '\\1\n', str(title), 0, re.DOTALL)
     t = re.sub(r'[(\']', ' ', t)
     t = re.sub(r'[\',)]', ' ', t)
            
     pp_filenodot= file_name.replace(".", "")

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
       #ax.xaxis.set_ticks_position('none')
       ax.yaxis.set_ticks_position('none')
             
      # Change the labels to the off-black
     ax.xaxis.label.set_color('#262626')
     ax.yaxis.label.set_color('#262626')

     if not os.path.exists(save_dir): os.makedirs(save_dir)

     plt.savefig('%s/EMBRACE_%s_%s_%s_%s_notitle.png' % (save_dir, area, pp_filenodot, ls, type_of_plot),
                 format='png', bbox_inches='tight')
                  
     fig.autofmt_xdate()
     

     #plt.title('\n'.join(wrap('%s' % (t.title()), 1000,replace_whitespace=False)), fontsize=16)
      #plt.show()
  
     #plt.savefig('%s/EMBRACE_Diurnal__monsoon_trough_%s_%s_%s.png' % (save_dir, pp_filenodot, ls, type_of_plot),
              #   format='png', bbox_inches='tight')
     #plt.close()

     #except Exception, e:
      #print e


if __name__ == '__main__':
   main()
       



