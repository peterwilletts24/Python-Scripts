import matplotlib
from matplotlib import rc
from matplotlib.font_manager import FontProperties
from matplotlib import rcParams
from matplotlib import cm

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import cm as cm_base

import matplotlib.pyplot as plt
import matplotlib.colors as colors

from matplotlib.patches import Polygon

from matplotlib.colors import from_levels_and_colors

rc('font', family = 'serif', serif = 'cmr10')
rc('text', usetex=True)

rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True
rcParams['font.family']='serif'
rcParams['font.serif']='cmr10'
rcParams['font.weight']='normal'
rcParams['text.color']='#262626'

import numpy as np

divisor=10  # for lat/lon rounding
parallels = np.arange(0.,90,divisor)
meridians = np.arange(0.,360., divisor)

ticks= np.arange(int(min_contour),int(max_contour)+tick_gap,tick_gap)

clevs = np.linspace(min_contour, max_contour,256)
midpoint=0
midp = np.mean(np.c_[clevs[:-1], clevs[1:]], axis=1)

vals = np.interp(midp, [min_contour, midpoint, max_contour], [0, 0.5, 1])
cols = plt.cm.RdBu_r(vals)

clevs_extend = np.linspace(min_contour, max_contour,254)
cmap, norm = from_levels_and_colors(clevs_extend, cols, extend='both')
