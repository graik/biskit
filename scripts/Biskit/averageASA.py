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
##

##
## last $Author$
## $Date$
## $Revision$
##

from Biskit import molUtils as MOU
import Biskit.mathUtils as MAU
from Biskit import tools as T
import os, os.path
import Numeric as N
from Biskit.PDBDope import PDBDope
from Biskit.PDBModel import PDBModel
import time, sys
import glob

def _use():
    print """
  averageASA.py - a script that collects the average (of 500 structures)
                  molecular surface (MS) and solvent accessible surface
                  (AS) for all 20 amino acids in a GLY-XXX-GLY tripeptide.

 Syntax:  AverageASA.py -i |template.pdb| -r |path| [ -l |str| ]

 Options: -i     file, pdb tripeptide template file
          -r     path, calculation root folder (many directories with be
                    created in this folder)
          -c     str, label for the three result dictionaries
           
 Result:  4 dictionaries AS, AS_sd, MS and MS_sd, written to root folder

EXAMPLE of a template file:
file name: 1gxg_input.pdb
ATOM      1  N   GLY A   1     106.354   3.981  -4.419  1.00  0.00      A
ATOM      2  CA  GLY A   1     107.169   5.048  -5.007  1.00  0.00      A
ATOM      3  C   GLY A   1     107.793   5.907  -3.912  1.00  0.00      A
ATOM      4  O   GLY A   1     107.238   6.030  -2.828  1.00  0.00      A
ATOM      5  N   XXX A   2     108.915   6.522  -4.175  1.00  0.00      A
ATOM      6  CA  XXX A   2     109.581   7.404  -3.178  1.00  0.00      A
ATOM      7  C   XXX A   2     108.699   8.597  -2.819  1.00  0.00      A
ATOM      8  O   XXX A   2     107.921   9.074  -3.644  1.00  0.00      A
ATOM      9  N   GLY A   3     108.824   9.073  -1.587  1.00  0.00      A
ATOM     10  CA  GLY A   3     108.030  10.213  -1.132  1.00  0.00      A
ATOM     11  C   GLY A   3     108.919  11.272  -0.489  1.00  0.00      A
ATOM     12  OT1 GLY A   3     109.737  10.907   0.339  1.00  0.00      A
ATOM     13  OT2 GLY A   3     108.770  12.432  -0.835  1.00  0.00      A
END


Default options:
"""
    o = defaultOptions()
    for key, value in o.items():
        print "\t-",key, "\t",value
        
    sys.exit(0)


def defaultOptions():
    return {'i':'1gxg_input.pdb',
            'r':'~/GLY_pep/',
            'l':''}

            
def randomPeptides(template, base_folder):
    """
    genarates 500 peptide orientations of each amino acid
    """
    for k in MOU.aaAtoms.keys():
        dir = 'GLY-%s-GLY'%(k)

        ## create derectory and copy generic template
        os.makedirs(base_folder + dir)
        os.system('cp %s %s/.'%(base_folder + template, base_folder + dir) )

        ## enter directort and greate template
        os.chdir( base_folder + dir )
        pattern = 's/XXX/%s/g'%(k)
        os.system( 'perl -i.bak -p -e %s %s'%(pattern, template) )

        ## generate xplor script then build model
        os.system( 'pdb2xplor.py -i %s'%template )
        os.system( 'xplor < %s_input_generate.inp'%template[:4] )

        ## create a run pcr folder
        pcrFolder = base_folder + dir + '_pcr'
        os.makedirs( pcrFolder )

        ## run PCR-MD
        os.chdir( base_folder )
        os.system('runPcr.py -t %s -r %s_pcr -h localhost'%(dir, dir) )


        ## don't run too many calculations at once
        time.sleep(100) # let previous calculation start
        nr_calc = os.system("ps awx | egrep '%sGLY-.*GLY_pcr/pcr_00' | wc -l"\
                            %base_folder )
        while nr_calc >= 2:
            time.sleep(15)


def randomSurfaces( base_folder, label ):
    """
    calculate surfaces for all peptides and return the
    average and SD
    """
    ## container for results and standard deviations
    MS,    AS    = {}, {}
    MS_sd, AS_sd = {}, {}

    ## loop over peptide directories
    for k in MOU.aaAtoms.keys():
        dir = 'GLY-%s-GLY_pcr/pcr_00'%(k)
        fLst = glob.glob( dir + '/*.pdb')

        msLst = []
        asLst = []

        ## loop over pdb files for each peptide
        for f in fLst:

            ## load peptide and remove waters and hydrogens
            m = PDBModel( base_folder + f )
            m = m.compress( m.maskProtein() * m.maskHeavy() )

            ## add surface data
            try:
                d = PDBDope( m )
                d.addSurfaceRacer( probe=1.4 )

                ## remove tailing GLY
                m = m.compress( m.res2atomMask([0,1,0]) )

                ## collect surface data for each peptide
                msLst += [ m.profile('MS') ]
                asLst += [ m.profile('AS') ]
                       
            except:
                print 'Failed calculating exposure for GLY-%s-GLY'%(k)
                print '\t and file %s'%f
                

        ## get result dictionary for peptide
        msDic = {}
        asDic = {}
        msDic_sd = {}
        asDic_sd = {}

        j = 0
        atoms =  [ a['name'] for a in m.atoms ]
        for n in atoms:
            msDic[n]    = N.average(msLst)[j]
            asDic[n]    = N.average(asLst)[j]
            msDic_sd[n] = MAU.SD( msLst )[j]
            asDic_sd[n] = MAU.SD( asLst )[j]
            j += 1

        MS[ k ] = msDic
        AS[ k ] = asDic
        MS_sd[ k ] = msDic_sd
        AS_sd[ k ] = asDic_sd

    return MS, AS, MS_sd, AS_sd

###########################
# MAIN
###########################

if len(sys.argv) < 3:
    _use()

options = T.cmdDict( defaultOptions() )

## where to run the calculation
base_folder = T.absfile( options['r'] )+'/'

## template gly-xxx-gly
template = options['i']

## label the result dictionaty files
label = '_' + options['l']

## create random prptides from template
randomPeptides( template, base_folder )

## collect average surfaces
MS, AS, MS_sd, AS_sd = randomSurfaces( base_folder, label )

## save dictionary with all 20 amino acids
T.Dump( MS, base_folder + 'MS%s.dic'%label)
T.Dump( AS, base_folder + 'AS%s.dic'%label )
T.Dump( MS_sd, base_folder + 'MS_sd%s.dic'%label )
T.Dump( AS_sd, base_folder + 'AS_sd%s.dic'%label )

