"""

Load npy xy, plot and save


"""

import os, sys

import matplotlib

matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!

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

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 14}

matplotlib.rc('font', **font)

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

top_dir='/nfs/a90/eepdw/Data/Rain_Land_Sea_Diurnal'

pp_file = 'avg.5216'

lon_max = 101.866 
lon_min = 64.115
lat_max= 33.
lat_min=-6.79

trmm_dir = '/nfs/a90/eepdw/Data/Observations/Satellite/TRMM/Diurnal/'
trmm_file = "trmm_diurnal_average_lat_%s_%s_lon_%s_%s.npz" % (lat_min,lat_max, lon_min, lon_max)

#############

# Make own time x-axis
utc_to_local=datetime.timedelta(hours=5, minutes=30)
d = matplotlib.dates.drange(datetime.datetime(2011, 8, 21, 6,30)+utc_to_local, datetime.datetime(2011, 8, 22, 6, 30)+utc_to_local, timedelta(hours=1))

formatter = matplotlib.dates.DateFormatter('%H:%M')

def main():
 #experiment_ids = ['djznw', 'djzny', 'djznq', 'djzns', 'dkjxq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq', 'dkbhu', 'djznu', 'dkhgu' ] # All 12
 experiment_ids_p = ['djznw', 'djzny', 'djznq', 'dklzq', 'dkmbq', 'dkjxq' ] # Most of Params
 experiment_ids_e = ['dklwu', 'dklyu', 'djzns', 'dkbhu', 'djznu', 'dkhgu'] # Most of Explicit
 #experiment_ids_p = ['dkmbq', 'dklzq' ] 
 #experiment_ids_e = ['dklwu', 'dklyu', 'djznu']

#experiment_ids = ['djzny', 'djznq', 'djzns', 'djznw', 'dkjxq', 'dklyu', 'dkmbq', 'dklwu', 'dklzq' ] 

 #plt.ion()

 NUM_COLOURS = 15
 cmap=cm.get_cmap(cm.Set1, NUM_COLOURS)
 #cgen = (cmap(1.*i/NUM_COLORS) for i in range(NUM_COLORS))

 for ls in ['land', 'total']:
 #for ls in ['sea']:
  fig = plt.figure(figsize=(16,8))
  ax = fig.add_subplot(111)
  legendEntries=[]
  legendtext=[]

  plot_trmm = np.load('%s%s_%s' % (trmm_dir, ls, trmm_file))
 
  dates_trmm=[]
  p=[]
  for dp in plot_trmm['hour']:
      print dp
      if ((int(dp)<23) & (int(dp)>=6)):
          dates_trmm.append(datetime.datetime(2011, 8, 21, int(dp), 0))
          p.append(plot_trmm['mean'][plot_trmm['hour']==dp])
      if ((int(dp)>=0) & (int(dp)<=6)):
          dates_trmm.append(datetime.datetime(2011, 8, 22, int(dp), 0))
          p.append(plot_trmm['mean'][plot_trmm['hour']==dp])
              
  #print dates_trmm
  a = np.argsort(dates_trmm,axis=0)

  d_trmm = np.array(dates_trmm)[a]
  pl = (np.array(p)[a])
  #pl=np.sort(pl,axis=1)
  
  l, = plt.plot_date(d_trmm+utc_to_local, pl, label='TRMM', linewidth=2*1.5, linestyle='-', marker='', markersize=2, fmt='', color='#262626')

  legendEntries.append(l)
  legendtext.append('TRMM')

  #land
  l0=plt.legend(legendEntries, legendtext,title='', frameon=False, loc=9, bbox_to_anchor=(0.31, 0,1, 1))
  #sea
  #l0=plt.legend(legendEntries, legendtext,title='', frameon=False, loc=9, bbox_to_anchor=(0, 0,1, 1))
  # Change the legend label colors to almost black
  texts = l0.texts
  for t in texts:
    t.set_color('#262626')

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
          colour = cmap(1.*2/NUM_COLOURS)
          linewidth=0.5
          linestylez='--'
   if ((experiment_id=='djznq') or (experiment_id=='dkjxq')):
          print experiment_id
          colour = cmap(1.*6/NUM_COLOURS)
          linewidth=0.8
          if (experiment_id=='djznq'):
              linestylez='--'
          if (experiment_id=='dkjxq'):
              linestylez=':'
              
   if ((experiment_id=='dklzq') or (experiment_id=='dklwu')):
          print experiment_id
          colour = cmap(1.*8/NUM_COLOURS)
          linewidth=1
          if (experiment_id=='dklzq'):
              linestylez='--'
          if (experiment_id=='dklwu'):
              linestylez='-'
   if ((experiment_id=='dklyu') or (experiment_id=='dkmbq')):
          print experiment_id
          colour = cmap(1.*10/NUM_COLOURS)
          linewidth=1.3
          if (experiment_id=='dkmbq'):
              linestylez='--'
          if (experiment_id=='dklyu'):
              linestylez='-'

   try:
      plotnp = np.load('%s/%s/%s/%s_%s_rainfall_diurnal_np_domain_constrain.npy' % (top_dir, expmin1, experiment_id, pp_file, ls))

      if (ls != 'total'):
          l, = plt.plot_date(d, plotnp[0]*3600, label=model_name_convert_legend.main(experiment_id), linewidth=linewidth*2.5, linestyle=linestylez, marker='', markersize=2, fmt='', color=colour)
      else:
          l, = plt.plot_date(d, plotnp*3600, label=model_name_convert_legend.main(experiment_id), linewidth=linewidth*2.5, linestyle=linestylez, marker='', markersize=2, fmt='', color=colour) 
      legendEntries.append(l)
      legendtext.append('%s' % (model_name_convert_legend.main(experiment_id)))
  

   except Exception, e:
      print e
      pass

  #Land
  l1=plt.legend(legendEntries, legendtext, title='Parametrised', loc=9, frameon=False, bbox_to_anchor=(0, 0,1, 1))
  #Sea
  #l1=plt.legend(legendEntries, legendtext, title='Parametrised', loc=9, frameon=False, bbox_to_anchor=(-0.255, 0,1, 1))
  
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


      
      if (ls != 'total'):
          l, = plt.plot_date(d, plotnp[0]*3600, label='%s' % (model_name_convert_legend.main(experiment_id)), linewidth=linewidth*1.5, linestyle=linestylez, marker='', markersize=2, fmt='', color=colour)
          
      else:
          l, = plt.plot_date(d, plotnp*3600, label='%s' % (model_name_convert_legend.main(experiment_id)), linewidth=linewidth*1.5, linestyle=linestylez, marker='', markersize=2, fmt='', color=colour)
      legendEntries.append(l)
      legendtext.append('%s' % (model_name_convert_legend.main(experiment_id)))

   except Exception, e:
      print e
      pass

  #Land
  l2=plt.legend(legendEntries, legendtext, title='Explicit', loc=9, frameon=False, bbox_to_anchor=(0.155, 0,1, 1))
  #Sea
  #l2=plt.legend(legendEntries, legendtext, title='Explicit', loc=9, frameon=False, bbox_to_anchor=(-0.1, 0,1, 1))
  
  plt.gca().add_artist(l1)
  plt.gca().add_artist(l0)
  plt.gca().xaxis.set_major_formatter(formatter)

  # Change the legend label colors to almost black
  texts = l2.texts
  for t in texts:
    t.set_color('#262626')

  plt.xlabel('Time (local)')
  plt.ylabel('mm/h')

  title="Domain Averaged Rainfall - %s" % ls

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

  datemin = d.min()
  datemax = d.max()
  ax.set_xlim(datemin, datemax)


  # Change the labels to the off-black
  ax.xaxis.label.set_color('#262626')
  ax.yaxis.label.set_color('#262626')


  if not os.path.exists('/nfs/a90/eepdw/Figures/EMBRACE/Diurnal/'): os.makedirs('/nfs/a90/eepdw/Figures/EMBRACE/Diurnal/')
  plt.savefig('/nfs/a90/eepdw/Figures/EMBRACE/Diurnal/%s_%s_latlon_dkbhu_notitle_big_font.png' % (pp_filenodot, ls), format='png', bbox_inches='tight')

  plt.title('\n'.join(wrap('%s' % (t.title()), 1000,replace_whitespace=False)), fontsize=16, color='#262626')
  #plt.show()
  
  plt.savefig('/nfs/a90/eepdw/Figures/EMBRACE/Diurnal/%s_%s_latlon_dkbhu_big_font.png' % (pp_filenodot, ls), format='png', bbox_inches='tight')
  plt.close()


if __name__ == '__main__':
   main()

   
       



