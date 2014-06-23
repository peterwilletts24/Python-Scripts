#!/usr/bin/env python
#   Copyright (C) 2004 CCLRC & NERC( Natural Environment Research Council ).
#   This software may be distributed under the terms of the
#   Q Public License, version 1.0 or later. http://ndg.nerc.ac.uk/public_docs/QPublic_license.txt

"""
na2cdms.py
==========

Module holds the na2cdms function to convert a NASA Ames file to a CDMS
file (NetCDF).

"""

helpMessage="""

na2nc
=====

Converts a NASA Ames file to a NetCDF file.

Usage
=====

   na2nc.py [<options>] -i <infilename> -o <outfilename> 

Where
-----

    infilename  - name of input file (NASA Ames).
    outfilename - name of output file (NetCDF).
    options     - list of options:
                    -t <time_units_string>  (such as "hours since 2003-04-30 10:00:00")
    
"""

# Imports from python standard library
import sys
sys.path.append(r"..")
import os
import re
import time

# Import from nappy package
import general
import cdmsMap

def na2cdms(nafile, ncfile, time_units=None, time_warning="yes", rename_variables={}, rules=None):
    print "Reading data from: %s" % nafile
    file=general.openNAFile(nafile)
    if file.FFI in (2110, 2160, 2310): 
	print """\nERROR: Cannot convert NASA Ames File Format Index (FFI) %s to NetCDF. 
No mapping implemented yet.""" % file.FFI
        return 0

    print "Writing output NetCDF file: %s\n" % ncfile
    file.toCdmsFile(ncfile, time_units=time_units, time_warning=time_warning, rename_variables=rename_variables)
    
    print "\nNetCDF file written successfully: %s" % ncfile

    return 1

if __name__=="__main__":

    args=sys.argv[1:]
    if len(args)<4:
        print helpMessage
        print "Incorrect number of arguments used."
        sys.exit()
	
    time_units=None
    time_warning="yes"
    rename_variables={}
    
    for arg in args:
        if arg=="-i":
	    infile=args[args.index(arg)+1]
	elif arg=="-o":
	    outfile=args[args.index(arg)+1]
	elif arg=="-t":
	    time_units=args[args.index(arg)+1]
	elif arg=="-n":
	    time_warning="no"
	elif arg=="-r":
	    renamer=args[args.index(arg)+1].split("=")	    
            rename_variables[renamer[0]]=renamer[1]

    na2cdms(infile, outfile, time_units=time_units, time_warning=time_warning, rename_variables=rename_variables) 

