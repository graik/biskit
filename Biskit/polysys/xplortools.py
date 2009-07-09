
from Biskit import PDBModel
import psfGen
from restools import linearize
import numpy as N
from vectors import norm

def linkProteins(a = None,residues_a = [], b = None,residues_b = [],linkseq = "",distance = 0):
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
    
    
    linearize(a, a.xyz[-1]-cm1,[0,0,1])
    linearize(b, b.xyz[0]-cm2,[0,0,-1])
    
    a.writePdb("lol1.pdb")
    b.writePdb("lol2.pdb")
    
    ## Translation
    ## - Real distance is from a's ending residue to b's starting residue
    distance = distance+norm(cm1)+norm(cm2)
    b.xyz += N.array([0,0,distance]) 
    
    ## Link creation
    link = PDBModel()

    for k in PDBModel.PDB_KEYS:
        link[k] = 5 * ['']  

    link['name'] = 5 * ['CA']
    link['residue_name'] = 5 * ['ALA']
    link['residue_number'] = N.arange( 5 ) + 1
    link['serial'] = N.arange( 5 ) + 1

    link.xyz = N.zeros( (5,3), float )
    
    ## Concatenation
    
    ## - Residue renumbering
    a.renumberResidues() ## a starts from 1
    last = a.atoms['residue_number'][-1]
    
    first = last+1
    
    b.renumberResidues(start= first)
    
    ## - Change chains and segments
    
    
    final = a.concat(b)
    final.writePdb("lol.pdb")


linkProteins(a = PDBModel('2AWT.pdb') , b= PDBModel('2CQJ.pdb'),distance = 70 )