
from Biskit import PDBModel
import psfGen
from restools import orientVectors
import numpy as N
from vectors import norm,normalized
from wormlikechain import WormLikeChainModel
from Biskit.molUtils import single2longAA

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
    
    
    orientVectors(a, a.xyz[-1]-cm1,[0,0,1])
    orientVectors(b, b.xyz[0]-cm2,[direction[0]*-1,direction[1]*-1,direction[2]*-1])
    
    ## Go to the origin!
    a.xyz -= cm1
    b.xyz -= cm2
    
    a.writePdb("lol1.pdb")
    b.writePdb("lol2.pdb")
    
    ## Translation
    ## - Real distance is from a's ending residue to b's starting residue
    total_distance = distance+norm(cm1)+norm(cm2)
    b.xyz += N.array([direction[0]*total_distance,direction[1]*total_distance,direction[2]*total_distance]) 
    
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
    
    ## Make triangular Calpha-trace
    linkdist_2 = linkdist / 2.
    
    dist_2 = distance / 2.
    
    
    
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
    
    a.atoms['segment_id'] = ["A"]*len(a.atoms['chain_id'])
    link.atoms['segment_id'] = ["L"]*len(link.atoms['chain_id'])
    b.atoms['segment_id'] = ["B"]*len(b.atoms['chain_id'])
    
    final = a.concat(link).concat(b)
    final.report()
    final.writePdb("lol.pdb")


linkProteins(a = PDBModel('2AWT.pdb') , b= PDBModel('2CQJ.pdb'),distance = 70 )