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
##
##
##
## load model dictionary and report
## * average pairwise rmsd
## * average rmsd to free structure (assumed to be model 1)
## * average rmsd to bound structure
## model 1 is not included in the calculations

## to run for all complexes:
## for c in c??; do
##   echo $c
##   cd $c/dock_pcr_multi_0919
##   a_model_rms.py pcr_rec/????_models.dic pcr_lig/????_models.dic ../com_wet/ref.complex
##   cd ~/interfaces
## done

from Biskit import *
from Biskit.tools import *
from Biskit.mathUtils import *
from Numeric import *


def interface_mask( refCom ):
    """
    -> interface mask for rec and lig of given complex
    """
    cont = refCom.resContacts( cutoff=4.5, force=1 )

    ## interface mask for rec and lig
    if_lig = refCom.lig_model.res2atomMask( N.sum( cont, 0) )
    if_lig = if_lig * refCom.lig_model.maskHeavy()

    if_rec = refCom.rec_model.res2atomMask( N.sum( cont, 1) )
    if_rec = if_rec * refCom.rec_model.maskHeavy()

    return if_rec, if_lig


def cast( drec, dlig, refCom ):
    """
    equalize atom content of set of free rec and lig with reference complex
    drec - { int : PDBModel } model dictionary
    dlig -       ~
    refCom - Complex, reference
    """
    rec_i_m, rec_i_c = drec[1].compareAtoms( refCom.rec_model )
    lig_i_m, lig_i_c = dlig[1].compareAtoms( refCom.lig_model )

    refCom.rec_model.keep( rec_i_c )
    refCom.lig_model.keep( lig_i_c )
    refCom.contacts = None

    for m in drec.values(): m.keep( rec_i_m )
    for m in dlig.values(): m.keep( lig_i_m )


def rms_pw_and_to_free( d, mask=None, bb=0, sc=0 ):
    """
    Get pw rms of models in fdic and their average rms to the free structure
    (which is assumed to be in position 1)
    d    - { int : PDBModel }
    mask - array of 0||1, atom mask (e.g. for interface)
    -> {..}
    """
    models = [ d[k] for k in range( 2, 12 ) ]

    t = Trajectory( models )
    if mask is None: mask = t.ref.maskHeavy()

    ## only consider backbone or side chains
    if bb:
        mask = mask * t.ref.maskBB()
    if sc:
        mask = mask * logical_not( t.ref.maskBB() )

    ## average pairwise RMS of all 10 snapshots
    pw = aboveDiagonal( t.pairwiseRmsd(mask) )
    pw_avg, pw_sd = average( pw )[0], SD( pw )[0]

    ## average RMSD to free (reference) structure
    t.fit( mask=mask, ref=d[1], prof='rms_ref' )
    rms_to_free = average( t.profile('rms_ref') )
    rms_to_free_sd = SD(t.profile('rms_ref') )

    return {'pdb':models[0].pdbCode,
            'pw' :pw_avg, 'pw_sd':pw_sd,
            'rms_f':rms_to_free, 'rms_f_sd':rms_to_free_sd }


def rms_to_bound( drec, dlig, ref, bb=0, sc=0 ):
    """
    Get average interface rms of rec and lig to reference complex
    Models in rec and lig and ref must be atom-casted to each other.
    drec, dlig   - { int : PDBModel }
    ref     - Complex
    """
    if_rec, if_lig = interface_mask( ref )

    trec = Trajectory( [ drec[k] for k in range( 2, 12 ) ] )
    tlig = Trajectory( [ dlig[k] for k in range( 2, 12 ) ] )

    if bb:
        if_rec = if_rec * trec.ref.maskBB()
        if_lig = if_lig * tlig.ref.maskBB()
    if sc:
        if_rec = if_rec * logical_not( trec.ref.maskBB() )
        if_lig = if_lig * logical_not( tlig.ref.maskBB() )

    trec.fit( mask=if_rec, ref=ref.rec(), prof='rms_if' )
    tlig.fit( mask=if_lig, ref=ref.lig(), prof='rms_if' )

    rms_rec = average( trec.profile( 'rms_if' ) )
    rms_rec_sd = SD( trec.profile('rms_if' ) )

    rms_lig = average( tlig.profile( 'rms_if' ))
    rms_lig_sd = SD( tlig.profile('rms_if'))

    return ( {'rms_b' : rms_rec, 'rms_b_sd' : rms_rec_sd },
             {'rms_b' : rms_lig, 'rms_b_sd' : rms_lig_sd } )


def rms_to_bound_2( drec, dlig ):
    """
    use saved values
    """
    recs = [ drec[k] for k in range( 2, 12 ) ]
    ligs = [ dlig[k] for k in range( 2, 12 ) ]

    rms_rec = average( [ m.info['rmsIF_to_bound'] for m in recs ] )
    rms_rec_sd = SD( [ m.info['rmsIF_to_bound'] for m in recs ] )

    rms_lig = average( [ m.info['rmsIF_to_bound'] for m in ligs ] )
    rms_lig_sd = SD( [ m.info['rmsIF_to_bound'] for m in ligs ] )

    return ( {'rms_b' : rms_rec, 'rms_b_sd' : rms_rec_sd },
             {'rms_b' : rms_lig, 'rms_b_sd' : rms_lig_sd } )
    

def report( d ):
    """
    d - { 'pdb':.., 'rms_f': .. }
    """
    print '%(pdb)s %(pw)5.1f %(pw_sd)4.2f  %(rms_f)5.1f %(rms_f_sd)4.2f  %(rms_b)5.1f %(rms_b_sd)4.2f' % d
    

###### MAIN #########

try:
    fdic_rec = sys.argv[1] ## model dict
    fdic_lig = sys.argv[2] 
    frefCom  = sys.argv[3] ## ref complex
except:
    print "\na_model_rms.py rec_models.dic lig_models.dic ref.complex\n"

errWriteln('Loading...')
d_rec = load( fdic_rec )
d_lig = load( fdic_lig )
refCom= load( frefCom ) 

errWrite('pairwise...')
## all heavy atom rms spread
result_rec = rms_pw_and_to_free( d_rec )
result_lig = rms_pw_and_to_free( d_lig )

## ## DEBUG ##

## rms_rec = [ rms_to_bound_2( rec, d_lig[1], refCom )[0] for rec in d_rec.values() ]

###########

cast( d_rec, d_lig, refCom )

## interface only rms spread
## if_mask_rec, if_mask_lig = interface_mask( refCom )
## result_rec = rms_pw_and_to_free( d_rec, if_mask_rec )
## result_lig = rms_pw_and_to_free( d_lig, if_mask_lig )

errWrite('\nto bound...')
to_bnd_rec, to_bnd_lig = rms_to_bound( d_rec, d_lig, refCom )

result_rec.update( to_bnd_rec )
result_lig.update( to_bnd_lig )

result_rec['pdb'] = d_rec[1].pdbCode
result_lig['pdb'] = d_lig[1].pdbCode

report( result_rec )

report( result_lig )
