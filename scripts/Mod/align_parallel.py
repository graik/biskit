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

from Biskit.Mod.AlignerMaster import AlignerMaster
from Biskit.Mod.ValidationSetup import ValidationSetup as VS
from Biskit.Mod.TemplateSearcher import TemplateSearcher as TS
from Biskit.Mod.Aligner import Aligner as A
import Biskit.tools as T
import Biskit.PVM.hosts as hosts
import sys, os
import os.path as osp
import glob

def _use( o ):

    print """
Built multiple alignment for each project given directory (parallelised).
If run from within a standardized modeling/validation folder structure,
i.e from the project root where the folders templates, sequences, and
validation reside all options will be set by the script.
        
Syntax: align_parallel.py -d |list of folders| -h |hosts|
                         [-pdb |pdbFolder| -ft |fastaTemplates|
                          -fs |fastaSequences| -fta |fastaTarget|
                          -fe |ferror|]

Note:  pvm must be running on the local machine!

Options:
    -d    [str], list of project directory (full path)
    -h    int, number of hosts to be used
    -a    first add hosts to pvm
    -pdb  str, pdbFolder for the pdb *.alpha
    -ft   str, path to 'templates.fasta'
    -fs   str, path to 'nr.fasta'
    -fta  str, path to 'target.fasta'
    -fe   str, path to the error file for the AlignerMaster

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
    ## does current folder look like a main project folder?
    if osp.exists( f + TS.F_RESULT_FOLDER ):
        d += [f]

    ## check if t_coffee folders alredy exist
    r = []
    for i in d:
        if osp.exists( i + A.F_RESULT_FOLDER ):
            print 'T-Coffee output folder alredy exists in %s'%i
            s = raw_input('Overwrite folder? (y/N)')
            if not s: s='N'
            if not ( s[0] =='y' or s[0]=='Y' ):
                r += [i]
    ## have too remove objects backvards not to change indexes
    r.reverse()
    for j in r:
        d.remove(j)
        
    if len(d)==0:
        print 'Nothing to align. Exiting.'
        sys.exit(0)

    return d

if __name__ == '__main__':

    ## look for default cross-validation projects
    d = []

    options = T.cmdDict({'h':10, 'd':d})

    d=options['d'] or collect_project_folders( )
    
    if (options['d'] is None) or ('help' in options or '?' in options):
        _use( options )

    folders = T.toList(options['d'])
    hostNumber = int(options['h'])
    pdbFolder = options.get('pdb', None)
    fastaTemplates = options.get('ft', None)
    fastaSequences = options.get('fs', None)
    fastaTarget = options.get('fta', None)
    ferror = options.get('fe',None)
    show_output = 'w' in options

    print "Initialize Job queue.."

    master = AlignerMaster(hosts=hosts.cpus_all[ : hostNumber ],
                           folders=folders,
                           pdbFolder = pdbFolder,
                           fastaTemplates=fastaTemplates,
                           fastaSequences=fastaSequences,
                           fastaTarget=fastaTarget,
                           ferror=ferror,
                           show_output=show_output)
    
    print "Starting alignment jobs .."
    master.calculateResult( )

    ## check for completed jobs
    j=0
    for i in d:
        if osp.exists( i + A.F_FINAL +'.pir_aln' ):
            j+=1
        else:
            print 'ERROR:Alignment error in %s see Aligner.log for more info'%i
        if j == len(d):
            print '\nAll %i alignments completed sucessfully.'%j
