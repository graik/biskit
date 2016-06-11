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
## 

from Numeric import *
import sys
import os, re
from Biskit.tools import *
from string import *
import biggles
from Biskit.Dock.Complex import Complex as ProteinComplex
from Biskit.mathUtils import *
import pipes
from time import time

def _use():
    print """
    CAD (contact area difference) calculation by icmbrowser.
    The calculation is performed only for residues that are in contact 
      in the reference complex.
      
        cl  - complexList, has to contain info dictionary data for key
        ref - reference complex

Default values:
    """
    default = _defOptions()
    for k in default:
        print "\t",k,"\t",default[k]
        
    sys.exit(0)


def _defOptions():
    return { 'cl':'complexes_cont.cl',
             'ref':'../../ref.complex'}


def nameFromPath( path ):
    """
    Extract a nice name from the absolute complexes_grouped.cg path
    """
    path = string.split( path, '/' )
    name = ''
    for p in path:
        if re.match( '^[a,c][0-9]{2}', p ):
            name += re.findall( '^[a,c][0-9]{2}', p )[0] + '-'
        if re.match(  '^dock_[0-9,a-z,_]{4}', p ):
            try:
                name += re.findall( '^dock_[0-9,a-z,_]{4}.+', p )[0]
            except:
                name += re.findall( '^dock_[0-9,a-z,_]{4}', p )[0]
            
    return name

    
def get_models( clst , ref_com ):
    rec_models = {}
    lig_models = {}
    for c in clst:
        if not c.info['model1'] in rec_models.keys():
            rec_models[c.info['model1']]= c.rec()
        if not c.info['model2'] in lig_models.keys():
            lig_models[c.info['model2']]= c.lig()

    rec_model_f_ref = rec_models[1]
    lig_model_f_ref = lig_models[1]

    rec_model_b_ref = ref_com.rec()
    lig_model_b_ref = ref_com.lig_model

    ## make sure ref and models have same atoms and atom order
    ## this should always be the true
    rec_m_ind, rec_f_ind = rec_models[2].compareAtoms( rec_model_f_ref )
    lig_m_ind, lig_f_ind = lig_models[2].compareAtoms( lig_model_f_ref )

    if rec_m_ind != rec_f_ind:
        print 'WARNING! Free receptor reference (xray) and'
        print 'receptor models have different atom order!!'
        rec_model_f_ref = rec_model_f_ref.take( rec_f_ind )

        for k in rec_models.keys()[1:]:
            rec_models[k] = rec_models[k].take( rec_m_ind )
            lig_models[k] = lig_models[k].take( lig_m_ind )
            
    if lig_m_ind != lig_f_ind:
        print 'Free ligand reference (xray) and '
        print 'ligand models have different atom order!!'
        lig_model_f_ref = lig_model_f_ref.take( lig_f_ind )

        for k in rec_models.keys()[1:]:
            rec_models[k] = rec_models[k].take( rec_m_ind )
            lig_models[k] = lig_models[k].take( lig_m_ind )

    ## make sure models and bound have same atoms and atom order
    rec_b2f_ind, rec_f2b_ind = rec_model_b_ref.compareAtoms( rec_model_f_ref )
    lig_b2f_ind, lig_f2b_ind = lig_model_b_ref.compareAtoms( lig_model_f_ref )

    for k in rec_models.keys():
        rec_models[k] = rec_models[k].take( rec_f2b_ind )
        lig_models[k] = lig_models[k].take( lig_f2b_ind )

    rec_model_b_ref = rec_model_b_ref.take( rec_b2f_ind )
    lig_model_b_ref = lig_model_b_ref.take( lig_b2f_ind )
            
    return rec_models, lig_models, rec_model_b_ref, lig_model_b_ref


def pw_rmsd( models_dic, ref_b, contacts ):
    """
    models_dic - model dictionary, where the free xray structure has key 1
    ref_b      - model, bound reference
    contacts   - atom contact mask
    -> dictionary, {1:[rmsd_b, rmsd_f], 2: ... }
    """
    rmsd = {}
    for k in models_dic.keys():
        rmsd_b = ref_b.rms( models_dic[k], contacts )
        rmsd_f = models_dic[1].rms( models_dic[k], contacts  )
        rmsd[k] = [ rmsd_b, rmsd_f ]
    return rmsd


def runICM( inp_lst, icm_log ):
    """
    inp - str, content of inp file
    -> int, run time in s
    """
    start_time = time()

    p = pipes.Template()

    ## '--' -> use STDIN, STDOUT
    p.append( '/usr/bin/icmbrowser', '--' ) 

    ##open STDIN for writing, connect STDOUT to file
    stdin = p.open( icm_log, 'w' ) 

    for inp in inp_lst:
            stdin.write( inp )

    stdin.close()

    return time() - start_time


def zeroChainId( model ):
    """
    set chainId to ''
    """
    #for a in model.getAtoms():
        #a['chain_id']=''
    model['chain_id'] = len(model) * ['']
    
def parse_icm( filename ):
    """
    Extract CAD value from icm-log file.
    filname - str, path to icm log file
    """
    result = []
    lines = open( filename ).readlines()

    if  len(lines) == 0:
        print '\nERROR: ICM result file empty'
        return

    for i in range( len(lines) ):  
        if re.match('^icm/def> 1.8\*Cad.+', lines[i]):
            result += [ float( string.strip( lines[i+1] ) ) ]

    return result
    
  
def pw_cad(  models_dic, ref_b, contacts ):
    """
    models_dic - model dictionary, where the free xray structure has key 1
    ref_b      - model, bound reference
    contacts   - atom contact mask
    -> dictionary, {1:[rmsd_b, rmsd_f], 2: ... }
    """
    cad = {}

    # residue index of contacting residues
    resIdx = nonzero( ref_b.atom2resMask( contacts ) )
    resIdx = ''.join([ str(i+1)+',' for i in resIdx ])[:-1]

    # bound reference
    ref_b_file = tempfile.mktemp('icm_bound_pdb')
    ref_b.renumberResidues()
    #for a in ref_b.getAtoms():
        #a['chain_id']=''
    ref_b['chain_id'] = [''] * len(ref_b)
    ref_b.writePdb( ref_b_file )

    # free reference (model 1)
    ref_f_file = tempfile.mktemp('icm_free_pdb')
    ref_f = models_dic[1].clone()
    ref_f.renumberResidues()
    #for a in ref_f.getAtoms():
        #a['chain_id']=''
    ref_f['chain_id'] = [''] * len( ref_f )
    ref_f.writePdb( ref_f_file )
    
    for k in models_dic.keys():
        # log file on disk
        log_file = tempfile.mktemp('icm_log')    

        # write temp pdb model to disk
        m_file = tempfile.mktemp('icm_pdb')

        models_dic[k].renumberResidues()
        #for a in models_dic[k].getAtoms():
            #a['chain_id']=''
        models_dic[k]['chain_id'] = [''] * len( models_dic[k] )
        
        models_dic[k].writePdb( m_file )

        inp_lst= ['call _startup \n'] + \
                 ['read pdb unix cat %s \n'%ref_b_file] + \
                 ['copy a_ "bound" \n'] + \
                 ['read pdb unix cat %s \n'%ref_f_file] + \
                 ['copy a_ "free" \n'] + \
                 ['read pdb unix cat %s \n'%m_file] + \
                 ['1.8*Cad(a_1./%s a_3./%s) \n'%(resIdx, resIdx)] + \
                 ['1.8*Cad(a_2./%s a_3./%s) \n'%(resIdx, resIdx)]
        
        runICM( inp_lst, log_file )
        cad[k] = parse_icm( log_file )
    
        os.unlink( m_file )
        os.unlink( log_file )
        
    os.unlink( ref_b_file )
    os.unlink( ref_f_file )

    return cad


def test():
    options = _defOptions()
    dir = '/home/Bis/johan/interfaces/c05/dock_multi_0919/hex1008/'
    options['ref'] = '/home/Bis/johan/interfaces/c05/com_wet/ref.complex'
    options['cl'] = dir +'complexes_cont.cl'
    return options
  

###############
# main function

if len(sys.argv) < 2:
    _use()
    
options = cmdDict( _defOptions() )
#    options = test()

## paths ...
cFile = os.path.abspath( options['cl'] )
dock_dir = os.path.split( absfile(options['cl']))[0]

## file names for output ...
outDir = os.path.dirname(cFile) + '/model_vs_rmsd'
if not os.path.exists( outDir ):
    os.mkdir( outDir )
    
cName = nameFromPath( cFile )
outName = outDir + '/' + cName + '_model-CAD.txt'
outName_over = outDir + '/' + cName + '_model-CAD_overwiev.txt'

## load
flushPrint('Loading complex list \n')
cList = load( options['cl'] )
flushPrint('Loading bound complex \n')
ref_com = load( options['ref'] )

## get models
rec_models, lig_models, rec_b_ref, lig_b_ref = get_models( cList, ref_com )

## contacts
ref_com_b_cast = ProteinComplex( rec_b_ref, lig_b_ref )
ref_b_resCont = ref_com_b_cast.resContacts()

rec_aContacts = rec_b_ref.res2atomMask( sum( ref_b_resCont, 1 ) )
lig_aContacts = lig_b_ref.res2atomMask( sum( ref_b_resCont, 0 ) )

## get pw_rmsds
rec_rmsd_cont = pw_rmsd( rec_models, rec_b_ref, rec_aContacts )
lig_rmsd_cont = pw_rmsd( lig_models, lig_b_ref, lig_aContacts )

rec_rmsd_bb = pw_rmsd( rec_models, rec_b_ref, rec_b_ref.maskBB() )
lig_rmsd_bb = pw_rmsd( lig_models, lig_b_ref, lig_b_ref.maskBB() )

## get pw cad
rec_cad = pw_cad(  rec_models, rec_b_ref, rec_aContacts )
lig_cad = pw_cad(  lig_models, lig_b_ref, lig_aContacts )


## print data per docking
out_file = open( outName, 'w')
head_0 = '#\t' + cName
head_1 = '# %11s %30s'%('rec','lig')
head_2 = '#%12s %10s %8s %16s %10s %8s'%\
         ('rms_bb','rms_cont','CAD','rms_bb','rms_cont','CAD')

out_file.write( head_0 + '\n' )
out_file.write( head_1 + '\n' )
out_file.write( head_2 + '\n' )

for r in rec_models.keys():
    for l in lig_models.keys():
        rec_rms_b_bb = rec_rmsd_bb[r][0]
        rec_rms_f_bb = rec_rmsd_bb[r][1]
        rec_rms_b_cont = rec_rmsd_cont[r][0]
        rec_rms_f_cont = rec_rmsd_cont[r][1]
        rec_cad_b = rec_cad[r][0]
        rec_cad_f = rec_cad[r][1]
        lig_rms_b_bb = lig_rmsd_bb[l][0]
        lig_rms_f_bb = lig_rmsd_bb[l][1]
        lig_rms_b_cont = lig_rmsd_cont[l][0]
        lig_rms_f_cont = lig_rmsd_cont[l][1]
        lig_cad_b = lig_cad[l][0]
        lig_cad_f = lig_cad[l][1]

        line = '%4i %4.2f %4.2f %4.2f %4.2f %5.2f %5.2f %4i %4.2f %4.2f %4.2f %4.2f %5.2f %5.2f'\
               %(r, rec_rms_b_bb, rec_rms_f_bb,
                 rec_rms_b_cont, rec_rms_f_cont, rec_cad_b, rec_cad_f,
                 l, lig_rms_b_bb, lig_rms_f_bb,
                 lig_rms_b_cont, lig_rms_f_cont, lig_cad_b, lig_cad_f )
        out_file.write( line + '\n')

out_file.close()


## print overview
out_file = open( outName_over, 'w')
head_0 = '#\t' + cName
head_1 = '# %20s %30s'%('rec','lig')
head_2 = '#%12s %10s %8s %16s %10s %8s'%\
         ('rms_bb','rms_cont','CAD','rms_bb','rms_cont','CAD')

out_file.write( head_0 + '\n' )
out_file.write( head_1 + '\n' )
out_file.write( head_2 + '\n' )

for i in rec_models.keys():
    rec_rms_b_bb = rec_rmsd_bb[i][0]
    rec_rms_f_bb = rec_rmsd_bb[i][1]
    rec_rms_b_cont = rec_rmsd_cont[i][0]
    rec_rms_f_cont = rec_rmsd_cont[i][1]
    rec_cad_b = rec_cad[i][0]
    rec_cad_f = rec_cad[i][1]
    lig_rms_b_bb = lig_rmsd_bb[i][0]
    lig_rms_f_bb = lig_rmsd_bb[i][1]
    lig_rms_b_cont = lig_rmsd_cont[i][0]
    lig_rms_f_cont = lig_rmsd_cont[i][1]
    lig_cad_b = lig_cad[i][0]
    lig_cad_f = lig_cad[i][1]

    line = '%4i %4.2f %4.2f %4.2f %4.2f %5.2f %5.2f %4i %4.2f %4.2f %4.2f %4.2f %5.2f %5.2f'\
           %(i, rec_rms_b_bb, rec_rms_f_bb,
             rec_rms_b_cont, rec_rms_f_cont, rec_cad_b, rec_cad_f,
             i, lig_rms_b_bb, lig_rms_f_bb,
             lig_rms_b_cont, lig_rms_f_cont, lig_cad_b, lig_cad_f )
    out_file.write( line + '\n')

out_file.close()
