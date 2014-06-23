#########################
#
# Take ship observations times, latitudes and longitudes and plot as animation
#
# Created by Peter Willetts on 27/11/2013
#
###########################

import os, sys
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.animation as animation
from matplotlib import collections
from datetime import datetime,timedelta
import numpy as np
import cPickle as pickle

# Import data, sorted into domain area, and date, by ship I.D - from marine__observations_read.py


plt_nam, plt_lat, plt_lon, plt_tim = pickle.load(open('marine_obs_times_and_positions.p', 'rb')) 

r_length = len(plt_tim)
    # Set figure and metadata

    #parallels = np.arange(0.,90,10.)

    #meridians = np.arange(0.,360.,10.)

    #m = Basemap(projection='cyl',\
           # llcrnrlat=-10.,urcrnrlat=30.,\
            #llcrnrlon=60.,urcrnrlon=105.,\
            #rsphere=6371200.,resolution='l',area_thresh=10000)# draw coastlines, state and country boundaries, edge of map.
    #m.drawcoastlines()
    #m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    #m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10) 
    #self.parallels = np.arange(0.,90,10.)

    #self.meridians = np.arange(0.,360.,10.)
    #i=0
    #x, y = m(plt_lon,plt_lat)
    #scat = m.scatter(x[i],y[i],3,marker='o',color='red')
    #fig =  plt.figure(figsize=(8,8)) 
   
    #ax = fig.add_axes([0.1,0.1,0.8,0.8])


#FLOOR = -10
#CEILING = 10

class AnimatedScatter(object):
    def __init__(self):
        global r
        r=0
        self.stream = self.data_stream()

        self.fig = plt.figure()
        #self.fig.canvas.mpl_connect('draw_event',self.forceUpdate)
        #self.ax = self.fig.add_axes([0.1,0.1,0.8,0.8])
        self.fig, self.ax = plt.subplots()
        self.ani = animation.FuncAnimation(self.fig, self.update, interval=500, 
                                       init_func=self.setup_plot)
       
    #def change_angle(self):
     #   self.angle = (self.angle + 1)%360

    #def forceUpdate(self, event):
      #  self.scat.changed()

    def setup_plot(self):
        
        self.basemap = self.basemap()

        X = next(self.stream)
        #print X
        self.scat = self.m.scatter(X[0], X[1],marker='o',color='red')

        #self.ax.set_xlim3d(FLOOR, CEILING)
        #self.ax.set_ylim3d(FLOOR, CEILING)

        return self.scat,

    def data_stream(self):
        global r
        r=0
       
        data = self.m(plt_lon[0], plt_lat[0])
        plt.title(str(plt_tim[r]))
       
        #print data
        while True:
            
            r += 1
            #print r
            if r <= r_length :
             data = self.m(plt_lon[r], plt_lat[r])
             
            #print data
            yield data

    def update(self, i):
        #global r
        #r += 1
        #print r
        data = next(self.stream)
        self.scat.remove()
        self.scat = self.m.scatter(data[0], data[1],marker='o',color='red')
        plt.title(str(plt_tim[r]))
        return self.scat,

    def show(self):
        plt.show()

    def basemap(self):
        #global  m
        
        self.m = Basemap(projection='cyl', llcrnrlat=-10.,urcrnrlat=30., llcrnrlon=60.,urcrnrlon=105., rsphere=6371200.,resolution='l',area_thresh=10000)
        
        # draw coastlines, state and country boundaries, edge of map.
        
        self.m.drawcoastlines()
        self.parallels = np.arange(0.,90,10.)

        self.meridians = np.arange(0.,360.,10.)
   
        self.m.drawparallels(self.parallels,labels=[1,0,0,0],fontsize=10)
        self.m.drawmeridians(self.meridians,labels=[0,0,0,1],fontsize=10)

        return self.m,


if __name__ == '__main__':
   a = AnimatedScatter()
   a.ani.save("Marine_observations_date_time_location.mp4", writer='mencoder', fps=5, dpi=72)
   #a.show()
