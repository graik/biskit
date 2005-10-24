#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 2 of the
## License, or any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.
##
## You find a copy of the GNU General Public License in the file
## license.txt along with this program; if not, write to the Free
## Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

## Contributions: Olivier PERIN, Raik Gruenberg
## last $Author$
## last $Date$
## $Revision$

from Biskit.Mod.ModelMaster import ModelMaster
from Biskit.Mod.ValidationSetup import ValidationSetup
import Biskit.tools as T
import Biskit.hosts as hosts
import sys, os
import os.path as osp
import glob

def _use( o ):

    print """
Syntax: model_parallel.py -d |list of folders| -h |host|
                       [-fta |fastaTarget| -pir |f_pir|
                       -tf |template_folder| -sm |starting_model|
                       -em |ending_model| -fe |ferror|]

    pvm must be running on the local machine!

Result: Parallel modelling for each project directory given
        
Options:
        -d    [str], project directories  (default: ./validation/*)
        -h    int, number of hosts to be used  (default: 10)
        -fta  str, path to find 'target.fasta'
        -pir  str, alignment filename
        -tf   str, directories for input atom files
        -sm   int, index of the first model
        -em   int, index of the last model
        -fe   str, filename to output errors from the Slave

Default options:\
"""
    for key, value in o.items():
        print "\t-",key, "\t",value

    sys.exit(0)


if __name__ == '__main__':

    ## look for default cross-validation projects
    f = os.getcwd() + ValidationSetup.F_RESULT_FOLDER
    d = []
    if osp.exists( f ):
        d = glob.glob( f+'/*' )
    
    options = T.cmdDict({'h':10, 'd':d})

    if (options['d'] is None) or ('help' in options or '?' in options):
        _use( options )
                       
    folders = T.toList(options['d'])
    hostNumber = int(options['h'])

 
    fastaTarget = options.get('fta', None)
    
    f_pir = options.get('pir', None)
   
    template_folder = options.get('tf', None)

    starting_model = options.get('sm', None)
  
    ending_model = options.get('em', None)

    ferror = options.get('fe', None)


    print "Initialize Job queue.."

    master = ModelMaster(hosts=hosts.cpus_all[ : hostNumber ], folders=folders,
                           fastaTarget=fastaTarget, f_pir=f_pir, template_folder=template_folder,
                           starting_model=starting_model, ending_model=ending_model, ferror=ferror)
    master.calculateResult()

