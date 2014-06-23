#!/bin/env python
# mosaic geo-tiff Aster Global DEM
# create roi-pac and Gamma headers
# by Ran Novitsky Nof 2011

import Image,sys,os,glob,getopt,scipy.misc,time
from pylab import *

def usage():
  sys.exit("""
A script for reading Aster Global DEM files and converting to ROI-PAC dem
or Gamma DEM.

usage:

    """+ sys.argv[0] +""" [output] [-g] [-h]
    
    output        : DEM output file name  
                    no file extention needed.
    -g            : Produce Gamma DEM only (4 byte real - big endien format)
                    default is to produce ROI-PAC DEM only (2 byte integer
                    - small endien format)
    -h            : Print this help text.
    
    The code will look in the current directory for subdirectories containing 
    ASTER DEM files and will create a mosaicked DEM for use by ROI-PAC
    and Gamma interferometry software.
    
    ASTER data is available from:
    http://www.gdem.aster.ersdac.or.jp/
    
    Download and unzip the files in current directory.
    current directory should contain directories in the format:
    ASTGTM2_NYYEXXX5 
    where YY in latitude and XXX is longitude.
    Run the script.
    
    The script assumes that the files are adjuctant tiles. 
    Missing tiles will get zero values.
        
    Created by Ran Novitsky Nof, 2011
    last update: Mar 29, 2013
""")
ers2dtypes ={
"IEEE4ByteReal" : float32,
"IEEE8ByteReal" : float64,
"Signed32BitInteger" : int32,
"Signed16BitInteger" : int16,
"Signed8BitInteger" : int8,
"Unsigned32BitInteger" : uint32,
"Unsigned16BitInteger" : uint16,
"Unsigned8BitInteger" : uint8,
"Complex64"           : complex64
}

ByteOrders = {"LSBFirst" : "little","MSBFirst" : "big"}

class GDEM:
    """class for ASTER Global DEM
    name = '' 				# Name of file 
    path = ''				# The path to the file
    byteorder = "little" # byte order
    dim = 0.0e-00			# pixel dimetions	
    width = 0				# number of sampels in each row
    length = 0				# number of rows
    east = 0.				# longitude coordinates of upper left corner 
    north = 0.				# latitude coordinates of upper left corner 
    im = ''         # DEM TIFF image data
    data = array([])				# DEM data values holder 
    ."""
    __name__ = 'GDEM'
    def __array__(self):
      return self
    def __call__(self):
      return 'GDEM'
    def __init__(self,filename):
      self.path,self.name = os.path.split(filename)
      if self.path == '': self.path='.'  
      self.dim = 0.0e-00
      self.width = 0
      self.length = 0
      self.east = 0.
      self.north = 0.
      self.im = ''
      self.byteorder = sys.byteorder
      self.data = array([])
      self.getdata()
    def __str__(self):
      return 'name = '+repr(self.name)+'\npath = '+repr(self.path)+'\ndim = '+repr(self.dim)+'\nwidth = '+repr(self.width)+'\nlength = '+repr(self.length)+'\neast = '+repr(self.east)+'\nnorth = '+repr(self.north)+'\n'
    def getdata(self):
      if (os.path.exists(self.path+os.sep+self.name)):
        self.im = Image.open(self.path+os.sep+self.name)
        self.im.mode = "I"
        self.dim = self.im.tag.get(33550)[0]
        self.width = self.im.tag.get(256)[0]
        self.length = self.im.tag.get(257)[0]
        self.east = self.im.tag.get(33922)[3]
        self.north = self.im.tag.get(33922)[4]
        self.data = scipy.misc.fromimage(self.im)
    def ers_header(self):
      header='''DatasetHeader Begin 
        Version = "6.4"
        Name ="'''+self.name+'''.dem.ers"
        LastUpdated     = '''+time.strftime("%a %b %d %H:%M:%S GMT %Y",time.gmtime())+''' 
        DataSetType     = ERStorage 
        DataType        = Raster 
        ByteOrder       = '''+[b[0] for b in ByteOrders.items() if self.byteorder in b][0]+''' 
        CoordinateSpace Begin 
                Datum = "WGS84"
                Projection = "GEODETIC"
                CoordinateType  = EN 
                Rotation        = 0:0:0.0 
        CoordinateSpace End 
        RasterInfo Begin 
               CellType        = '''+[b[0] for b in ers2dtypes.items() if self.data.dtype.type in b][0]+''' 
               NullCellValue   = 0 
                CellInfo Begin 
                        Xdimension      = '''+str(self.dim)+''' 
                        Ydimension      = '''+str(self.dim)+'''
                CellInfo End 
                NrOfLines       = '''+str(int(self.length))+'''
                NrOfCellsPerLine        = '''+str(int(self.width))+''' 
                RegistrationCoord Begin 
                        Eastings        = '''+str(self.east-self.dim/2.)+''' 
                        Northings       = '''+str(self.north+self.dim/2.)+'''
                RegistrationCoord End 
               NrOfBands       =  1
               BandId Begin 
                        Value = "Elevation"
               BandId End
        RasterInfo End 
DatasetHeader End'''
      savetxt(self.path+os.sep+self.name+'.dem.ers',[header],'%s')        
    def roi_header(self):
      header='''WIDTH         '''+str(self.width)+'''
FILE_LENGTH   '''+str(self.length)+'''
XMIN          0
XMAX          '''+str(self.width-1)+'''
YMIN          0
YMAX          '''+str(self.length-1)+'''
X_FIRST       '''+str(self.east)+'''
Y_FIRST       '''+str(self.north)+'''
X_STEP        '''+str(self.dim)+'''
Y_STEP        '''+str(-self.dim)+'''
X_UNIT        degres
Y_UNIT        degres
Z_OFFSET      0
Z_SCALE       1
PROJECTION    LATLON
DATUM         WGS84
'''
      savetxt(self.path+os.sep+self.name+'.dem.rsc',[header],'%s')
    def gamma_header(self):
      header='''Gamma DIFF&GEO DEM/MAP parameter file
title: '''+self.name+'''
DEM_projection:     EQA
data_format:        REAL*4
DEM_hgt_offset:          0.00000
DEM_scale:               1.00000
width:                '''+str(self.width)+'''
nlines:               '''+str(self.length)+'''
corner_lat:     '''+str(self.north)+'''  decimal degrees
corner_lon:     '''+str(self.east)+'''  decimal degrees
post_lat:   -'''+str(self.dim)+'''  decimal degrees
post_lon:    '''+str(self.dim)+'''  decimal degrees

ellipsoid_name: WGS 84
ellipsoid_ra:        6378137.000   m
ellipsoid_reciprocal_flattening:  298.2572236

datum_name: WGS 1984
datum_shift_dx:              0.000   m
datum_shift_dy:              0.000   m
datum_shift_dz:              0.000   m
datum_scale_m:         0.00000e+00
datum_rotation_alpha:  0.00000e+00   arc-sec
datum_rotation_beta:   0.00000e+00   arc-sec
datum_rotation_gamma:  0.00000e+00   arc-sec
datum_country_list Global Definition, WGS84, World
'''
      savetxt(self.path+os.sep+self.name+'.dem.par',[header],'%s') 
    def save_GDEM(self):
      self.data.tofile(self.path+os.sep+self.name+'.dem')
      self.ers_header()
           
def get_dem_files():
  dem_files = glob.glob('AST*/AST*dem*')
  dem_files.sort()
  dems = []
  for dem in dem_files:
    dems = dems +[GDEM(dem)]
  return dems
  
def get_dems_extent(dems):
  return min([dem.east for dem in dems]),max([dem.east+dem.width*dem.dim for dem in dems]),min([dem.north-dem.length*dem.dim for dem in dems]),max([dem.north for dem in dems]),dems[0].dim

def stich(gdem,dem):
  gdem.data[int(round((gdem.north-dem.north)/dem.dim)):int(round((gdem.north-dem.north)/dem.dim+dem.length)),int(round((dem.east-gdem.east)/dem.dim)):int(round((dem.east-gdem.east)/dem.dim+dem.width))] = dem.data[:,:]

def mk_gdem(outfile):
  print "\tReading GDEM files..."
  dems = get_dem_files()
  gdem = GDEM(outfile)  
  west,east,south,north,dim = get_dems_extent(dems)
  gdem.dim = dim
  gdem.width = round((east-west)/dim)
  gdem.length = round((north-south)/dim)
  gdem.east = west
  gdem.north = north
  gdem.data = zeros((gdem.length,gdem.width),dtype=float32)
  print "\tMosaic GDEM files..."
  for dem in dems:
    stich(gdem,dem)
  gdem.im = scipy.misc.toimage(gdem.data)
  return gdem
  

if __name__=="__main__":
  print " ****** Mosaic ASTER Global DEM ******"
  try:
    opts,args = getopt.gnu_getopt(sys.argv[1:],'gh')
  except getopt.error, err:
    print str(err)
    usage()
  opts = dict(opts)
  if '-h' in opts: usage()
  try:
    outfile = args[0]
  except IndexError:
    print "Error: no output name"
    usage()
  gdem = mk_gdem(outfile) 
  if '-g' in opts:
    print "\tSaving to Gamma format..."
    gdem.data[where(gdem.data==0)]=0.0001
    gdem.data = gdem.data.byteswap().astype(float32)
    gdem.byteorder = 'big'
    gdem.gamma_header()
  else:
    print "\tSaving to ROI-PAC format..."
    gdem.data = gdem.data.astype(int16)
    gdem.roi_header()
  gdem.save_GDEM()  

  
  
  
