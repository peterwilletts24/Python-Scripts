# This serves as an intensive exercise of matplotlib's transforms
# and custom projection API. This example produces a so-called
# SkewT-logP diagram, which is a common plot in meteorology for
# displaying vertical profiles of temperature. As far as matplotlib is
# concerned, the complexity comes from having X and Y axes that are
# not orthogonal. This is handled by including a skew component to the
# basic Axes transforms. Additional complexity comes in handling the
# fact that the upper and lower X-axes have different data ranges, which
# necessitates a bunch of custom classes for ticks,spines, and the axis
# to handle this.

#import sys
#sys.path.append('/nfs/see-fs-01_users/eepdw/.local/lib/python2.7/site-packages/')

import numpy as np
from numpy import ma,array,linspace,log,cos,sin,pi,zeros,exp,arange,load,interp,pi,arctan2,sqrt,min,max,where
from matplotlib.axes import Axes
import matplotlib.transforms as transforms
import matplotlib.axis as maxis
import matplotlib.spines as mspines
import matplotlib.path as mpath
from matplotlib.projections import register_projection

import imp
imp.load_source('SoundingRoutines', '/nfs/see-fs-01_users/eepdw/python_scripts/Tephigram/Sounding_Routines.py')
from SoundingRoutines import *

#from skewt.thermodynamics import VirtualTemp,Latentc,SatVap,MixRatio,GammaW,\
#	VirtualTempFromMixR,MixR2VaporPress,DewPoint,Theta,TempK
#from skewt.thermodynamics import Rs_da, Cp_da, Epsilon
# The sole purpose of this class is to look at the upper, lower, or total
# interval as appropriate and see what parts of the tick to draw, if any.
class SkewXTick(maxis.XTick):
    def draw(self, renderer):
        if not self.get_visible(): return
        renderer.open_group(self.__name__)

        lower_interval = self.axes.xaxis.lower_interval
        upper_interval = self.axes.xaxis.upper_interval

        if self.gridOn and transforms.interval_contains(
                self.axes.xaxis.get_view_interval(), self.get_loc()):
            self.gridline.draw(renderer)

        if transforms.interval_contains(lower_interval, self.get_loc()):
            if self.tick1On:
                self.tick1line.draw(renderer)
            if self.label1On:
                self.label1.draw(renderer)

        if transforms.interval_contains(upper_interval, self.get_loc()):
            if self.tick2On:
                self.tick2line.draw(renderer)
            if self.label2On:
                self.label2.draw(renderer)

        renderer.close_group(self.__name__)


# This class exists to provide two separate sets of intervals to the tick,
# as well as create instances of the custom tick
class SkewXAxis(maxis.XAxis):
    def __init__(self, *args, **kwargs):
        maxis.XAxis.__init__(self, *args, **kwargs)
        self.upper_interval = 0.0, 1.0

    def _get_tick(self, major):
        return SkewXTick(self.axes, 0, '', major=major)

    @property
    def lower_interval(self):
        return self.axes.viewLim.intervalx

    def get_view_interval(self):
        return self.upper_interval[0], self.axes.viewLim.intervalx[1]


# This class exists to calculate the separate data range of the
# upper X-axis and draw the spine there. It also provides this range
# to the X-axis artist for ticking and gridlines
class SkewSpine(mspines.Spine):
    def _adjust_location(self):
        trans = self.axes.transDataToAxes.inverted()
        if self.spine_type == 'top':
            yloc = 1.0
        else:
            yloc = 0.0
        left = trans.transform_point((0.0, yloc))[0]
        right = trans.transform_point((1.0, yloc))[0]

        pts  = self._path.vertices
        pts[0, 0] = left
        pts[1, 0] = right
        self.axis.upper_interval = (left, right)


# This class handles registration of the skew-xaxes as a projection as well
# as setting up the appropriate transformations. It also overrides standard
# spines and axes instances as appropriate.
class SkewXAxes(Axes):
    # The projection must specify a name.  This will be used be the
    # user to select the projection, i.e. ``subplot(111,
    # projection='skewx')``.
    name = 'skewx'

    def _init_axis(self):
        #Taken from Axes and modified to use our modified X-axis
        self.xaxis = SkewXAxis(self)
        self.spines['top'].register_axis(self.xaxis)
        self.spines['bottom'].register_axis(self.xaxis)
        self.yaxis = maxis.YAxis(self)
        self.spines['left'].register_axis(self.yaxis)
        self.spines['right'].register_axis(self.yaxis)

    def _gen_axes_spines(self):
        spines = {'top':SkewSpine.linear_spine(self, 'top'),
                  'bottom':mspines.Spine.linear_spine(self, 'bottom'),
                  'left':mspines.Spine.linear_spine(self, 'left'),
                  'right':mspines.Spine.linear_spine(self, 'right')}
        return spines

    def _set_lim_and_transforms(self):
        """
        This is called once when the plot is created to set up all the
        transforms for the data, text and grids.
        """
        rot = 30

        #Get the standard transform setup from the Axes base class
        Axes._set_lim_and_transforms(self)

        # Need to put the skew in the middle, after the scale and limits,
        # but before the transAxes. This way, the skew is done in Axes
        # coordinates thus performing the transform around the proper origin
        # We keep the pre-transAxes transform around for other users, like the
        # spines for finding bounds
        self.transDataToAxes = self.transScale + (self.transLimits +
                transforms.Affine2D().skew_deg(rot, 0))

        # Create the full transform from Data to Pixels
        self.transData = self.transDataToAxes + self.transAxes

        # Blended transforms like this need to have the skewing applied using
        # both axes, in axes coords like before.
        self._xaxis_transform = (transforms.blended_transform_factory(
                    self.transScale + self.transLimits,
                    transforms.IdentityTransform()) +
                transforms.Affine2D().skew_deg(rot, 0)) + self.transAxes

    def add_dry_adiabats(self,T0,P,do_labels=True,**kwargs):
	# Added by Thomas Chubb
	P0=1000.
	T=array([ (st+273.15)*(P/P0)**(Rs_da/Cp_da)-273.15 for st in T0 ])
	labelt=[ (st+273.15)*1**(Rs_da/Cp_da) for st in T0 ]
	if kwargs.has_key('color'): 
	    col=kwargs['color']
	else: 
	    col='k'
	for tt,ll in zip(T,labelt):
	    self.plot(tt,P,**kwargs)
	    if do_labels:
		if (tt[8]>-50) and (tt[8]<20):
		    self.text(tt[8],P[8]+10,'%d'%(ll),fontsize=8,\
			    ha='center',va='bottom',rotation=-30,color=col,\
			    bbox={'facecolor':'w','edgecolor':'w', 'alpha':0.5})
	return T
    

    def add_moist_adiabats(self,T0,P,do_labels=True,**kwargs):
	# Added by Thomas Chubb
	T=array([lift_wet(st,P) for st in T0])
	if kwargs.has_key('color'): 
	    col=kwargs['color']
	else: 
	    col='k'
	for tt in T:
	    self.plot(tt,P,**kwargs)
	    # if (tt[-1]>-60) and (tt[-1]<-10):
	    if do_labels:
		self.text(tt[-1],P[-1],'%d'%tt[0],ha='center',va='bottom',\
			fontsize=8, bbox={'facecolor':'w','edgecolor':'w', 'alpha':0.5},\
                        color=col)

    def add_mixratio_isopleths(self,w,P,do_labels=True,**kwargs):
	# Added by Thomas Chubb
	e=array([P*ww/(.622+ww) for ww in w])
	T = 243.5/(17.67/log(e/6.112) - 1)
	if kwargs.has_key('color'): 
	    col=kwargs['color']
	else: 
	    col='k'

	for tt,mr in zip(T,w):
	    self.plot(tt,P.flatten(),**kwargs)
	    if do_labels:
		if (tt[0]>-45) and (tt[0]<15):
		    if mr*1000<1.:
			fmt="%4.1f"
		    else:
			fmt="%d"
		    self.text(tt[0],P[0],fmt%(mr*1000),\
			    color=col, fontsize=8,ha='center',va='bottom',\
			    bbox={'facecolor':'w','edgecolor':'w', 'alpha':0.5})

    def other_housekeeping(self,pmax, pmin, tmax, tmin, mixratio=array([])):
	# Added by Thomas Chubb
	self.yaxis.grid(True,ls='-',color='k',lw=0.5)
        
	# Plot x-grid lines instead of using xaxis.grid().
	# This is because xaxis.grid only plots skew lines 
	# that intersect the lower x axis... a possible fix
	# would be to use twinx() to plot upper and lower
	# grid lines but I think that's messy too.
	for TT in linspace(-100,100,21):
	    self.plot([TT,TT],[pmax,pmin],color='k',lw=0.5)
            if TT>=tmin and TT<=tmax:
                 self.text(TT,975,"%d" % TT,\
		           color='k', fontsize=12,ha='center',va='bottom',\
		           bbox={'facecolor':'w','edgecolor':'w', 'alpha':0.5})
	# self.set_ylabel('Pressure (hPa)')
	self.set_xlabel('Temperature (C)')

# Now register the projection with matplotlib so the user can select
# it.
register_projection(SkewXAxes)

def lift_parcel(self,startp,startt,startdp):
	"""Lift a parcel to discover certain properties.
	
	
	INPUTS:
	startp:  Pressure (hPa)
	startt:  Temperature (C)
	startdp: Dew Point Temperature (C)
	"""
	from numpy import interp

	assert startt>startdp,"Not a valid parcel. Check Td<Tc"
	Pres=linspace(startp,100,100)

	# Lift the dry parcel
	T_dry=(startt+273.15)*(Pres/startp)**(Rs_da/Cp_da)-273.15 

	# Mixing ratio isopleth
	starte=SatVap(startdp)
	startw=MixRatio(starte,startp*100)
	e=Pres*startw/(.622+startw)
	T_iso=243.5/(17.67/log(e/6.112)-1)

	# Solve for the intersection of these lines (LCL).
	# interp requires the x argument (argument 2)
	# to be ascending in order!
	P_lcl=interp(0,T_iso-T_dry,Pres)
	T_lcl=interp(P_lcl,Pres[::-1],T_dry[::-1])

	col=[.6,.6,.6]

	# zorder
	zo=4

	# Plot traces below LCL
	self.skewxaxis.plot(T_dry[Pres>=P_lcl],Pres[Pres>=P_lcl],color=col,lw=2,zorder=zo)
	self.skewxaxis.plot(T_iso[Pres>=P_lcl],Pres[Pres>=P_lcl],color=col,lw=2,zorder=zo)
	self.skewxaxis.plot(T_lcl,P_lcl,ls='',marker='o',mec=col,mfc=col,zorder=zo)

	# Now lift a wet parcel from the intersection point
	preswet=linspace(P_lcl,200)
	tempwet=lift_wet(T_lcl,preswet)

	# Plot trace above LCL
	self.skewxaxis.plot(tempwet,preswet,color=col,lw=2,zorder=zo)

	# Add text to sounding
	dtext ="Parcel:\n"
	dtext+="Ps:  %6.1fhPa\n"%startp
	dtext+="Ts:    %4.1fC\n"%startt
	dtext+="Ds:    %4.1fC\n"%startdp
	dtext+="Plcl:%6.1fhPa\n"%P_lcl
	dtext+="Tlcl:  %4.1fC"%T_lcl

	self.fig.text(0.1,0.895,dtext,fontname="monospace",va='top')

	return

def surface_parcel(self,mixdepth=125):
	"""Returns parameters for a parcel initialised by:
	1. Surface pressure (i.e. pressure of lowest level)
	2. Surface temperature determined from max(theta) of lowest <mixdepth> mbar
	3. Dew point temperature representative of lowest <mixdepth> mbar

	Inputs:
	mixdepth (mbar): depth to average mixing ratio over
	"""

	pres=self.data["pres"]
	temp=self.data["temp"]
	dwpt=self.data["dwpt"]

	# identify the layers for averaging
	layers=pres>pres[0]-mixdepth
	
	# parcel pressure is surface pressure
	pres_s=pres[0]

	# average theta over mixheight to give
	# parcel temperature
	thta_mix=Theta(temp[layers]+273.15,pres[layers]*100.).max()
	temp_s=TempK(thta_mix,pres_s*100)-273.15

	# average mixing ratio over mixheight
	vpres=SatVap(dwpt)
	mixr=MixRatio(vpres,pres*100)
	mixr_mix=mixr[layers].mean()
	vpres_s=MixR2VaporPress(mixr_mix,pres_s*100)

	# surface dew point temp
	dwpt_s=DewPoint(vpres_s)

	return pres_s,temp_s,dwpt_s

def lift_wet(startt,pres):
    #--------------------------------------------------------------------
    # Lift a parcel moist adiabatically from startp to endp.
    # Init temp is startt in C, pressure levels are in hPa    
    #--------------------------------------------------------------------

    temp=startt
    t_out=zeros(pres.shape);t_out[0]=startt
    for ii in range(pres.shape[0]-1):
	delp=pres[ii]-pres[ii+1]
 	temp=temp-100*delp*GammaW(temp+273.15,(pres[ii]-delp/2)*100)
	t_out[ii+1]=temp

    return t_out

def u_v_winds(wind_direction, wind_speed):
    wind_rad = np.radians(wind_direction)
    u_wind=-((wind_speed)*np.sin(wind_rad))
    v_wind=-((wind_speed)*np.cos(wind_rad))
    return u_wind,v_wind

if __name__ == '__main__':
    main()
