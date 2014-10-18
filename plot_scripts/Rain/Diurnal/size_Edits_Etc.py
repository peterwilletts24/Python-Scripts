font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 14}

matplotlib.rc('font', **font)

prop={'size':8}, #delete

# Make own time x-axis
utc_to_local=datetime.timedelta(hours=5, minutes=30)
d = matplotlib.dates.drange(datetime.datetime(2011, 8, 21, 6,30)+utc_to_local, datetime.datetime(2011, 8, 22, 6, 30)+utc_to_local, timedelta(hours=1))


#Sea

l0=plt.legend(legendEntries, legendtext,title='', frameon=False, loc=9, bbox_to_anchor=(0, 0,1, 1))

l1=plt.legend(legendEntries, legendtext, title='Parametrised', loc=9, frameon=False, bbox_to_anchor=(-0.255, 0,1, 1))

l2=plt.legend(legendEntries, legendtext, title='Explicit', loc=9, frameon=False, bbox_to_anchor=(-0.1, 0,1, 1))

#Land

l0=plt.legend(legendEntries, legendtext,title='', frameon=False, loc=9, bbox_to_anchor=(0.31, 0,1, 1))

l1=plt.legend(legendEntries, legendtext, title='Parametrised', loc=9, frameon=False, bbox_to_anchor=(0, 0,1, 1))

l2=plt.legend(legendEntries, legendtext, title='Explicit', loc=9, frameon=False, bbox_to_anchor=(0.155, 0,1, 1))

linewidth=linewidth*1.5


fig = plt.figure(figsize=(16,8))

datemin = d.min()
  datemax = d.max()
  ax.set_xlim(datemin, datemax)

  l, = plt.plot_date(d_trmm+utc_to_local, pl, label='TRMM', linewidth=2, linestyle='-', marker='', markersize=2, fmt='', color='#262626')

Update savefig name and filename

