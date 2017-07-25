#!/usr/bin/env python
## numpy-oldnumeric calls replaced by custom script; 09/06/2016
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
##
##

## Select frames from trajectory for docking


import Biskit.tools as T
import Biskit.mathUtils as MaU
from Biskit import TrajCluster, EnsembleTraj, PCRModel, molUtils
from Biskit.Dock import hexTools
from Biskit.EnsembleTraj import traj2ensemble
import Biskit.oldnumeric as N0

import os.path
import copy, sys

def _use():
    print """
selectModels: Select non-redundant frames from a trajectory dump them and put
them into a PDB file for HEX docking.

Syntax:  selectModels -i |traj.dat| -o |out_folder| [ -psf |psf_file|
                      -dic |out_model_dic|
                      -n |number| -ref
                      -co |TrajCluster.out| -a |atom1 atom2 atom..|
                      -s |startFrame| -e |endFrame| -step |frameSkip|
                      -id |chaiID]
                      -conv |convergence_diff| ]
                      
         i    - pickled Trajectory object
         dic  - alternative name for model.dic
         psf  - create PCRModels with psf file info
         ref  - add trajectory's reference model to dictionary and pdb if
                  a reference pdb file is given this will be used insted
         id   - set ligand and receptor chainID
         a    - atoms to use for clustering,
                default: C and roughly every second side chain heavy atom
         conv - float, convergence criterium [1e-11]

Result:  - n pickled PDBModels '|PDBCode|_|frame|.model' in out_folder
         - pickled TrajCluster if requested
         - |PDBCode|_model.dic with these n PDBModels indexed from 1 to n

Defaults:
"""
    default = defOptions()
    for k in default:
        print "\t",k,"\t",default[k]
        
    sys.exit(0)

def defOptions():
    return {'o':'.', 'n':'10',
            's':'0','step':'3',
            'conv':1e-11,
            'dic':None,
##            'hex':None,
##            'id':None
            }


def load( options ):

    f_traj = options['i']
    start = int( options['s'] )
    
    traj = T.load( f_traj )

    if traj.__class__ != EnsembleTraj:
        traj = traj2ensemble( traj )

    end = traj.lenFrames()
    if options.has_key('e'):
        end = int( options['e'] )

    step = int( options['step'] )

    traj = traj.takeFrames( range( start, end ) )
    if step != 1:
        traj = traj.thin( step )

    traj.fit()

    return TrajCluster( traj )

    
def selectedAtoms( model ):
    """
    -> atom mask with C, and roughly every second side chain heavy atom
    """
    red={'GLY':['C'],
            'ALA':['C','CB'],
            'VAL':['C','CB','CG2'],
            'LEU':['C','CB','CG'],
            'ILE':['C','CB','CG2'],
            'MET':['C','CB','SD'],
            'PRO':['C','CB','CD'],
            'PHE':['C','CB','CG','CZ'],
            'TRP':['C','CB','CD1','NE1','CE3','CZ3'],
            'SER':['C','OG'],
            'THR':['C','CB','CG2'],
            'ASN':['C','CB','OD1'],
            'GLN':['C','CB','CD','NE2'],
            'TYR':['C','CB','CG','CZ'],
            'CYS':['C','CB'],
            'LYS':['C','CB','CD','NZ'],
            'ARG':['C','CB','CD','CZ'],
            'HIS':['C','CB','ND1','CE1'],
            'ASP':['C','CG'],
            'GLU':['C','CB','CD']}

    result = []
    for res in model.resList():
        allowedAtoms = red[ res[0]['residue_name'] ]
        
        for a in res:
            result += [ (a['name'] in allowedAtoms) ]

    return result
    

def report( tc ):

    clTrajs = tc.memberTrajs()

    for i in range(0, tc.n_clusters ):
        
        t = clTrajs[i]
        rms = tc.avgRmsd( i, tc.aMask )

        names = [ '_'.join(T.stripFilename(s).split('_')[-2:])
                  for s in t.frameNames]

        print "%i <%4.2f +-%4.2f>: " % (i, rms[0],rms[1] ), names
        print

    tr = clTrajs[0].concat( *tuple( clTrajs[1:] ) )

    avgall = N0.average( MaU.aboveDiagonal( tr.pairwiseRmsd( tc.aMask ) ) )
    print "avg rms all: %4.2f" %  avgall


def dumpModel( m, options, fout ):

    if options.has_key('psf'):
        m = PCRModel( options['psf'], m )

    m.addChainFromSegid( verbose=0 ) # chain id removed by xplor
    m.removeRes( 'TIP3' )

    m.saveAs( fout )

    return m


def setChainID( m ):
    """
    set chainID for Hex pdb files
    """
    if options['id']:
        id = T.toList( options['id'] )
        cMap = m.chainMap()
        for chain in range( m.lenChains() ):
            idx = N0.nonzero( cMap == chain )
            for i in idx:
                m.atoms['chain_id'][i] = id[chain]


def rmsdLimitedClustering( tc, options, min_cluster=5, max_cluster=30,
                           rmsLimit=1.0 ):
    """
    Cluster iteratively until the average of all clusters meet
    the rmsd cutoff criteria.
    """

    allowedAtoms = T.toList( options.get('a',[]) )
    if allowedAtoms:
        mask = tc.traj.ref.mask( lambda a: a['name'] in allowedAtoms )
    else:
        mask = selectedAtoms( tc.traj.ref )


    n_cluster = tc.calcClusterNumber( min_clst=min_cluster,
                                      max_clst=max_cluster,
                                      rmsLimit=1.0 )

    options['n'] = n_cluster
    
    cluster( tc, options )
    


    
def cluster( tc, options ):
    
    n_cluster = int( options['n'] )

    allowedAtoms = T.toList( options.get('a',[]) )
    if allowedAtoms:
        mask = tc.traj.ref.mask( lambda a: a['name'] in allowedAtoms )
    else:
        mask = selectedAtoms( tc.traj.ref )
    
    saveIn = T.absfile( options['o'] ) + '/'
    conv = float( options['conv'] )
    
    tc.cluster( n_cluster, aMask=mask, converged=conv )

    ## collect center frame index for each cluster
    frames = [ members[0] for members in tc.memberFrames() ]

    result = tc.traj.takeFrames( frames ) ## trajectory of cluster centers

    model_dic = {}

    dic_index = 1

    if options.has_key('ref'):
        ## use user-provided reference structure
        if os.path.isfile( T.absfile(options['ref']) ):
            print '\nUsing user specified reference pdb'
            m = PDBModel( options['ref'] )
            m.remove( m.maskH2O() )
            
        ## use reference in trajectory
        else:
            print '\nUsing reference in trajectory' 
            m = tc.traj.ref

        m = dumpModel( m, options, saveIn+m.getPdbCode()+'_ref.model')
        ## add ref as first model in dictionary   
        model_dic[dic_index] = m
        dic_index += 1

    ## save the individual models and add them to the model dictionary
    for i in range(0, result.lenFrames() ):
        m = result.getPDBModel(i)

        m = dumpModel(m, options, saveIn +
                  T.stripFilename(result.frameNames[i]) +'.model' )

        model_dic[dic_index] = m
        dic_index += 1
        
        
    ## save model dictionary
    fdic = options['dic'] or m.getPdbCode() + '_models.dic'
    T.dump( model_dic, T.absfile( fdic ) )

## REDUNDANT CODE AS MULTIDOCK NOW WRITES THE HEX PDB FILES
##
##     ## save all models in the dictionary as HEX pdb files
##     for k in model_dic.keys():
##         m = model_dic[k]
        
##         ## remove hydrogens and sort atoms in standard order
##         m.remove( m.maskH() )
##         m = molUtils.sortAtomsOfModel(m)
##         setChainID( m )

##         ## save single hex pdbs
##         if options['hex']:
##             fhex = options['hex'] + '_%03d' %(k) 
##         else:
##             fhex = m.getPdbCode() + '_%03d_hex.pdb'%(k)
   
##         hexTools.createHexPdb_single( m, T.absfile( fhex ) )
       
#    fhex = options['hex'] or m.getPdbCode() + '_hex.pdb'
#    hexTools.createHexPdb( model_dic, T.absfile( fhex ) )
    
    return result


def test():
    options = defOptions()
    options['i'] = '/home/Bis/raik/data/tb/interfaces/c11/lig_pcr_00/traj.dat'
    options['step'] = '3'
    options['s'] = '0'
    options['o'] = '~johan/dock/scripts'
    options['ref'] = ''

    return options


if __name__ == '__main__':
    if len(sys.argv) < 2:
        _use()

##    options = test()
    options = T.cmdDict( defOptions() )
                       
    tc = load( options )

    r = cluster( tc, options )
##    r=rmsdLimitedClustering( tc, options )
    if options.has_key('co'):
        T.dump( tc, options['co'] )

    report( tc )
        
