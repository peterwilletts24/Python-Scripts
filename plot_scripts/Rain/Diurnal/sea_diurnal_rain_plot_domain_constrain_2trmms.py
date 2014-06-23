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

rc('text', usetex=True)

rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True

rc('font', family = 'serif', serif = 'cmr10')

import numpy as np

from datetime import timedelta
import datetime

import imp

import re
from textwrap import wrap

model_name_convert_legend = imp.load_source('util', '/nfs/see-fs-01_users/eepdw/python_scripts/modules/model_name_convert_legend.py')
#unrotate = imp.load_source('util', '/home/pwille/python_scripts/modules/unrotate_pole.py')

###############
# Things to change

pp_file = 'avg.5216'

trmm_dir = '/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/Diurnal/'

lon_max = 101.866 
lon_min = 64.115
lat_max= 33.
lat_min=-6.79

#lon_max = 116.
#lon_min = 30.5
#lat_max= 40.
#lat_min=-11.25

trmm_file = "trmm_diurnal_average_lat_%s_%s_lon_%s_%s.npz" % (lat_min,lat_max, lon_min, lon_max)
trmm_file2 = "trmm_diurnal_average_lat_%s_%s_lon_%s_%s_height_0.0_to_600.0.npz" % (lat_min,lat_max, lon_min, lon_max)
#############

# Make own time x-axis

d = matplotlib.dates.drange(datetime.datetime(2011, 8, 21, 6,30), datetime.datetime(2011, 8, 22, 6, 30), timedelta(hours=1))

formatter = matplotlib.dates.DateFormatter('%H:%M')

#print d

#times = matplotlib.dates.date2num(d)

top_dir='/nfs/a90/eepdw/Data/Rain_Land_Sea_Diurnal'

def main():
 #experiment_ids = ['djznw', 'djzny', 'djznq', 'djzns', 'dkjxq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq', 'dkbhu', 'djznu', 'dkhgu' ] # All 12
 experiment_ids_p = ['djznw', 'djzny', 'djznq', 'dklzq', 'dkmbq', 'dkjxq' ] # Most of Params
 experiment_ids_e = ['dklwu', 'dklyu', 'djzns', 'dkbhu', 'djznu', 'dkhgu'] # Most of Explicit
#experiment_ids = ['djzny', 'djznq', 'djzns', 'djznw', 'dkjxq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq' ] 

 #plt.ion()

 NUM_COLOURS = 15
 cmap=cm.get_cmap(cm.Set1, NUM_COLOURS)
 #cgen = (cmap(1.*i/NUM_COLORS) for i in range(NUM_COLORS))

 for ls in ['land']:
  plt.figure(figsize=(12,6))

  legendEntries=[]
  legendtext=[]

  plot_trmm = np.load('%s%s_%s' % (trmm_dir, ls, trmm_file))
  plot_trmm2 = np.load('%s%s_%s' % (trmm_dir, ls, trmm_file2))
  dates_trmm=[]
  p=[]

  dates_trmm2=[]
  p2=[]
  for dp in plot_trmm['hour']:
      print dp
      if ((int(dp)<23) & (int(dp)>=6)):
          dates_trmm.append(datetime.datetime(2011, 8, 21, int(dp), 0))
          p.append(plot_trmm['mean'][plot_trmm['hour']==dp])
      if ((int(dp)>=0) & (int(dp)<=6)):
          dates_trmm.append(datetime.datetime(2011, 8, 22, int(dp), 0))
          p.append(plot_trmm['mean'][plot_trmm['hour']==dp])
  for dp in plot_trmm2['hour']:
      print dp
      if ((int(dp)<23) & (int(dp)>=6)):
          dates_trmm2.append(datetime.datetime(2011, 8, 21, int(dp), 0))
          p2.append(plot_trmm2['mean'][plot_trmm2['hour']==dp])
      if ((int(dp)>=0) & (int(dp)<=6)):
          dates_trmm2.append(datetime.datetime(2011, 8, 22, int(dp), 0))
          p2.append(plot_trmm2['mean'][plot_trmm2['hour']==dp])              
  #print dates_trmm
  a = np.argsort(dates_trmm,axis=0)
  a2 = np.argsort(dates_trmm2,axis=0)
  
  d_trmm = np.array(dates_trmm)[a]
  d_trmm2 = np.array(dates_trmm2)[a2]

  pl = (np.array(p)[a])/3600
  pl2 = (np.array(p2)[a2])/3600

#pl=np.sort(pl,axis=1)
  
  l, = plt.plot_date(d_trmm, pl, label='TRMM', linewidth=2, linestyle='-', marker='', markersize=2, fmt='', color='black')

  legendEntries.append(l)
  legendtext.append('TRMM')

  l, = plt.plot_date(d_trmm, pl2, label='TRMM600m', linewidth=2, linestyle='-', marker='', markersize=2, fmt='', color='grey')

  legendEntries.append(l)
  legendtext.append('TRMM less than 600m')
  l0=plt.legend(legendEntries, legendtext,title='', frameon=False, prop={'size':8}, loc=9, bbox_to_anchor=(0.21, 0,1, 1))

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
          colour = cmap(1.*5/NUM_COLOURS)
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
      plotnp = np.load('%s/%s/%s/%s_%s_rainfall_diurnal_np_domain_constrain.npy' % (top_dir, expmin1, experiment_id, pp_file, ls))

      l, = plt.plot_date(d, plotnp[0], label=model_name_convert_legend.main(experiment_id), linewidth=linewidth, linestyle=linestylez, marker='', markersize=2, fmt='', color=colour)

      legendEntries.append(l)
      legendtext.append('%s' % (model_name_convert_legend.main(experiment_id)))
  
  
   except Exception, e:
      print e
      pass

  l1=plt.legend(legendEntries, legendtext, title='Parametrised', loc=9, frameon=False, prop={'size':8}, bbox_to_anchor=(0, 0,1, 1))

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
          colour = cmap(1.*5/NUM_COLOURS)
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
      plotnp = np.load('%s/%s/%s/%s_%s_rainfall_diurnal_np_domain_constrain.npy' % (top_dir, expmin1, experiment_id, pp_file, ls))


      #plotnp = np.sort(pnp, axis=1)
      l, = plt.plot_date(d, plotnp[0], label='%s' % (model_name_convert_legend.main(experiment_id)), linewidth=linewidth, linestyle=linestylez, marker='', markersize=2, fmt='', color=colour)

      legendEntries.append(l)
      legendtext.append('%s' % (model_name_convert_legend.main(experiment_id)))

   except Exception, e:
      print e
      pass
  l2=plt.legend(legendEntries, legendtext, title='Explicit', loc=9, frameon=False, bbox_to_anchor=(0.11, 0,1, 1), prop={'size':8})
  plt.gca().add_artist(l1)
  plt.gca().add_artist(l0)
  plt.gca().xaxis.set_major_formatter(formatter)

  plt.xlabel('Time (UTC)')
  plt.ylabel('kg/m^2/s')

  title="Domain Averaged Rainfall - %s" % ls

  t=re.sub('(.{68} )', '\\1\n', str(title), 0, re.DOTALL)
  t = re.sub(r'[(\']', ' ', t)
  t = re.sub(r'[\',)]', ' ', t)
  
  
  plt.title('\n'.join(wrap('%s' % (t.title()), 1000,replace_whitespace=False)), fontsize=16)
  plt.show()

  if not os.path.exists('/nfs/a90/eepdw/Rainfall_plots/Figures/Diurnal/'): os.makedirs('/nfs/a90/eepdw/Rainfall_plots/Figures/Diurnal/')
  #plt.savefig('/nfs/a90/eepdw/Rainfall_plots/Figures/Diurnal/%s_%s.png' % (pp_file, ls), format='png', bbox_inches='tight')
  plt.close()


if __name__ == '__main__':
   main()

   
       



