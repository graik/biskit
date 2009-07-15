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

def cleanPdb(m1):
    print m1
    ## cleanup structures
    m1 = m1.compress( m1.maskProtein() * m1.maskHeavy() )
    
    m1.remove( m1.indicesFrom( 'name', ['OXT'] ) )
    m1.remove( m1.indicesFrom( 'name', ['OT1'] ) )
    m1["name"][ m1.indicesFrom( 'name', ['OT2'] )] = 'O'
    
    m1 = PDBCleaner( m1 ).process()
    
    return m1
    
def linkProteins(a = None, b = None,link_seq = "GSGSGSGSA",direction = [0,0,1],distance = 0):
    """
    Links two proteins with a flexible peptide of repeating pattern linkseq.
    If distance is 0 then it tries to use the whole peptide.
    
    @param a: First model to be linked.
    @type a: PDBModel
    @param b: Second model to be linked.
    @type b: PDBModel
    @param residues_a: Residues to be used to get the center of mass of the first 
            operator.
    @type residues_a:list{int}
    @param residues_b: Residues to be used to get the center of mass of the second
            operator.
    @type residues_b: list{int}
    
    TODO
    """
    
    direction = N.array(normalized(direction))
    
    cm1 = a.centerOfMass()
    
    cm2 = b.centerOfMass()
    
    a = cleanPdb(a)
    b = cleanPdb(b)
    
    
    ## Go to the origin!
    a.xyz -= cm1
    b.xyz -= cm2
    
    
    ## Orient
    b_start = b.xyz[N.where(b.maskFrom( 'name', ['N'] ))[0][0]]
    a_ending = a.xyz[N.where(a.maskFrom( 'name', ['C'] ))[0][-1]]
    
    orientVectors(a, a_ending-cm1,direction)
    orientVectors(b, b_start-cm2,[direction[0]*-1,direction[1]*-1,direction[2]*-1])
    
    
    ## Link creation
    ## - Get contour from distance
    
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
    linkseq = ""
    for i in range(nres):
        linkseq += link_seq[i%len(link_seq)]
        
    print "Link sequence:", linkseq

    ## Concatenation
    ## - Change chains and segments
    a.atoms['chain_id'] = ["A"]*len(a.atoms['chain_id'])
    b.atoms['chain_id'] = ["A"]*len(b.atoms['chain_id'])
    
    a.atoms['segment_id'] = ["p1"]*len(a.atoms['chain_id'])
    b.atoms['segment_id'] = ["p1"]*len(b.atoms['chain_id'])
    
    a['serial'] = range( len(a) )
    b['serial'] = range( len(b) )
    
    open("01_sequence.txt","w").write(a.sequence()+linkseq+b.sequence())
    open("01_link_sequence.txt","w").write(linkseq)
    #~ open('01_resid_a.txt',"w").write("resid "+str(a.atoms['residue_number'][0])+":"+str(a.atoms['residue_number'][-1]))
    #~ open('01_resid_b.txt',"w").write("resid "+str(b.atoms['residue_number'][0])+":"+str(b.atoms['residue_number'][-1]))
    #~ open('01_resid_link.txt',"w").write("resid "+str(link.atoms['residue_number'][0])+":"+str(link.atoms['residue_number'][-1]))
    
    ## Create and 'link' the peptide
    os.system("xplor -py xplor_build_link.py")
    link = PDBModel("01_extended_link.pdb")
    link.xyz -= link.centerOfMass()
    link_end = link.xyz[N.where( link.maskFrom( 'name', ['OT1'] ))[0][0]] 
    link_start = link.xyz[N.where( link.maskFrom( 'name', ['HT2'] ))[0][0]] 
    orientVectors(link,link_end-link_start,direction)
    link_end = link.xyz[N.where( link.maskFrom( 'name', ['OT1'] ))[0][0]] 
    link_start = link.xyz[N.where( link.maskFrom( 'name', ['HT2'] ))[0][0]] 
    cleanPdb(link)
        
    ## - Residue renumbering
    last = a.atoms['residue_number'][-1] 
    first = last+ 1 + nres
    a.renumberResidues() ## a starts from 1
    link.renumberResidues(start = last+1)
    b.renumberResidues(start= first)
    
    ## Linking...
    a_last_C = a.xyz[N.where(a.maskFrom( 'name', ['C'] ))[0][-1]]
    b_first_N = b.xyz[N.where(b.maskFrom( 'name', ['N'] ))[0][0]]
    
    a.xyz += link_start - a_last_C
    b.xyz += link_end - b_first_N
    
    link.atoms['chain_id'] = ["A"]*len(link.atoms['chain_id'])
    link.atoms['segment_id'] = ["L"]*len(link.atoms['chain_id'])
    
    final = a.concat(link).concat(b)
    final.report()
    final =  cleanPdb(final)
    final.writePdb("final.pdb") 
    
    
linkProteins(a = PDBModel('2AWT.pdb') , b= PDBModel('2CQJ.pdb'),distance = 70 )

#~ link = PDBModel()
#~ link.xyz = N.zeros( (nres,3), float )

#~ for k in PDBModel.PDB_KEYS:
    #~ link[k] = nres * [PDBModel.PDB_KEYS[k][1]]  

#~ link['name'] = nres * ['CA']

#~ linkseq3 = single2longAA(linkseq)

#~ reslist = []
#~ for i in range(len(link['name'])):
    #~ reslist.append(linkseq3[i%len(linkseq)])
#~ link['residue_name'] = reslist
#~ link['residue_number'] = N.arange( nres ) + 1
#~ link['serial'] = N.arange( nres ) + 1
    