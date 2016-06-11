#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2016 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
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

from Biskit.Mod.ModelMaster import ModelMaster
from Biskit.Mod.ValidationSetup import ValidationSetup as VS
from Biskit.Mod.TemplateSearcher import TemplateSearcher as TS
from Biskit.Mod.Modeller import Modeller as M
from Biskit.Mod.TemplateFilter import TemplateFilter as TF
import Biskit.tools as T
import Biskit.PVM.hosts as hosts
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


def collect_project_folders( ):

    d=[]
    f = os.getcwd()
    if osp.exists( f + VS.F_RESULT_FOLDER ):
        d = glob.glob( f + VS.F_RESULT_FOLDER + '/????' )
    ## does current look like a main project folder?
    if osp.exists( f + TS.F_RESULT_FOLDER ):
        d += [f]

    ## check if a modeller folders alredy exist
    r = []  
    for i in d:
        if osp.exists( i + M.F_RESULT_FOLDER ):
            print 'Modeller output folder alredy exists in %s'%i
            s = raw_input('Overwrite folder? (y/N)')
            if not s: s='N'
            if not ( s[0] =='y' or s[0]=='Y' ):
                r += [i]
    ## have too remove objects backvards not to change indexes
    r.reverse()
    for j in r:
        d.remove(j)
        
    if len(d)==0:
        print 'Nothing to model. Exiting.'
        sys.exit(0)

    return d


if __name__ == '__main__':

    ## look for default cross-validation projects
    options = T.cmdDict({'h':10, 'd':[], 'zfilter':TF.Z_CUTOFF,'idfilter':TF.ID_CUTOFF})
    d=options['d'] or collect_project_folders( )

    if (options['d'] is None) or ('help' in options or '?' in options):
        _use( options )
                       
    folders = T.toList(options['d'])
    hostNumber = int(options['h'])
    fastaTarget = options.get('fta', None)
    f_pir = options.get('pir', None)
    template_folder = options.get('tf', None)
    starting_model = int(options.get('sm', 1))
    ending_model = int(options.get('em', 10))
    ferror = options.get('fe', None)
    windows = "w" in options

    print "Initialize Job queue.."

    master = ModelMaster(hosts=hosts.cpus_all[ : hostNumber ],
                         folders=folders,
                         fastaTarget=fastaTarget,
                         f_pir=f_pir,
                         template_folder=template_folder,
                         starting_model=starting_model,
                         ending_model=ending_model,
                         ferror=ferror,show_output=windows)
    
    master.calculateResult()

