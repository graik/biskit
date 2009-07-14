from restools import orientVectors
import numpy as N
from vectors import norm,normalized
from wormlikechain import WormLikeChainModel
from Biskit.molUtils import single2longAA
from Biskit.PDBModel import PDBModel 
import os


import Biskit.tools as T
import Biskit.molUtils as MU
from Biskit import PDBCleaner

def cleanPdbs(m1,m2):
    print m1,m2
    ## cleanup structures
    m1 = m1.compress( m1.maskProtein() * m1.maskHeavy() )
    m2 = m2.compress( m2.maskProtein() * m2.maskHeavy() )

    m1.remove( m1.indicesFrom( 'name', ['OXT'] ) )
    m2.remove( m2.indicesFrom( 'name', ['OXT'] ) )

    m1 = PDBCleaner( m1 ).process()
    m2 = PDBCleaner( m2 ).process()

    return m1,m2
    
def linkProteins(a = None,residues_a = [], b = None,residues_b = [],linkseq = "GSGSGSGSA",direction = [0,0,1],distance = 0):
    """
    Links two proteins with a flexible peptide of repeating pattern linkseq.
    If distance is 0 then it tries to use the whole peptide.
    
    @param a: First model to be linked.
    @type a: PDBModel
    @param b: First model to be linked.
    @type b: PDBModel
    @param residues_a: Residues to be used to get the center of mass of the first 
            operator.
    @type residues_a:list{int}
    @param residues_b: Residues to be used to get the center of mass of the second
            operator.
    @type residues_b: list{int}
    
    TODO
    """
    if(residues_a != []):
        cm1 = a.takeResidues(residues_a).centerOfMass()
    else:
        cm1 = a.centerOfMass()
        
    if(residues_b != []):
        cm2 = b.takeResidues(residues_b).centerOfMass()
    else:
        cm2 = b.centerOfMass()
    
    a,b = cleanPdbs(a,b)
    
    
    a_ending = a.resModels()[-1].centerOfMass()
    b_ending = b.resModels()[-1].centerOfMass()
    
    
    ## Orient
    orientVectors(a, a_ending-cm1,direction)
    orientVectors(b, b_ending-cm2,[direction[0]*-1,direction[1]*-1,direction[2]*-1])
    
    a_ending = a.resModels()[-1].centerOfMass()
    b_ending = b.resModels()[-1].centerOfMass()
        
    ## Go to the origin!
    a.xyz -= cm1
    b.xyz -= cm2
    
    a.writePdb("lol1.pdb")
    b.writePdb("lol2.pdb")
    
    
    ## Link creation
    ## - Get distance
    
    wlc = WormLikeChainModel()
    p = 4.; x = distance 
        
    wlc.p = p;wlc.x=x;
    wlc.T = 298.;wlc.F = 10e-12

    linkdist = wlc.getContourLengthFromX()
    
    print linkdist
    
    linkdist = 40 ## in A 
    
    ## - Get number of residues (3.5A for N to N)
    nres = int( linkdist / 3.5 )
    print nres

    ## - Create link structure
    link = PDBModel()
    link.xyz = N.zeros( (nres,3), float )
    
    for k in PDBModel.PDB_KEYS:
        link[k] = nres * [PDBModel.PDB_KEYS[k][1]]  
    
    link['name'] = nres * ['CA']
    
    linkseq3 = single2longAA(linkseq)

    reslist = []
    for i in range(len(link['name'])):
        reslist.append(linkseq3[i%len(linkseq)])
    link['residue_name'] = reslist
    link['residue_number'] = N.arange( nres ) + 1
    link['serial'] = N.arange( nres ) + 1
    
        

    ## Concatenation
    
    ## - Residue renumbering
    last = a.atoms['residue_number'][-1] 
    first = last+ 1 + nres
    a.renumberResidues() ## a starts from 1
    link.renumberResidues(start = last+1)
    b.renumberResidues(start= first)
    
    ## - Change chains and segments
    a.atoms['chain_id'] = [" "]*len(a.atoms['chain_id'])
    link.atoms['chain_id'] = [" "]*len(link.atoms['chain_id'])
    b.atoms['chain_id'] = [" "]*len(b.atoms['chain_id'])
    
    a.atoms['segment_id'] = ["p1"]*len(a.atoms['chain_id'])
    link.atoms['segment_id'] = ["L"]*len(link.atoms['chain_id'])
    b.atoms['segment_id'] = ["p1"]*len(b.atoms['chain_id'])
    
    a['serial'] = range( len(a) )
    b['serial'] = range( len(b) )
    
    final = a.concat(b)
    final.report()
    final.writePdb("aux_in.pdb")
    
    a.writePdb( '01_domain1.pdb' )
    b.writePdb( '01_domain2.pdb' )
    
    open("01_sequence.txt","w").write(a.sequence()+link.sequence()+b.sequence())
    open("01_link_sequence.txt","w").write(link.sequence())
    open('01_resid_a.txt',"w").write("resid "+str(a.atoms['residue_number'][0])+":"+str(a.atoms['residue_number'][-1]))
    open('01_resid_b.txt',"w").write("resid "+str(b.atoms['residue_number'][0])+":"+str(b.atoms['residue_number'][-1]))
    open('01_resid_link.txt',"w").write("resid "+str(link.atoms['residue_number'][0])+":"+str(link.atoms['residue_number'][-1]))
    
    os.system("xplor -py xplor_rebuild_script.py")
    #~ link = PDBModel("aux_link.pdb")
    #~ ap = link.centerOfMass()
    #~ bp = link.resModels()[-1].centerOfMass()
    #~ orientVectors(link,bp-ap,direction)
    
    
linkProteins(a = PDBModel('2AWT.pdb') , b= PDBModel('2CQJ.pdb'),distance = 70 )