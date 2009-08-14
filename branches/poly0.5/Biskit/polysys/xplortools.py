from restools import orientVectors, orientPlanes
import numpy as N
from vectors import norm,normalized
from wormlikechain import WormLikeChainModel
from Biskit.molUtils import single2longAA
from Biskit.PDBModel import PDBModel 
import os
from xtender import Xtender 

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
    orientVectors(b, b_start-cm2,direction*-1)
    
    
    ## Link creation
    ## - Get contour from distance
    
    wlc = WormLikeChainModel()
    p = 4.; x = distance 
        
    wlc.p = p;wlc.x=x;
    wlc.T = 298.;wlc.F = 10e-12

    linkdist = wlc.getContourLengthFromX()
    
    print linkdist
    
    linkdist = 700 ## in A 
    
    ## - Get number of residues (3.5A for N to N)
    nres = int( linkdist / 3.5 )
    print nres

    ## - Create link structure
    
    linkseq = ""
    for i in range(nres):
        linkseq += link_seq[i%len(link_seq)]
        
    #~ open("01_sequence.txt","w").write(a.sequence()+linkseq+b.sequence())
    #~ open("01_link_sequence.txt","w").write(linkseq)
    #~ os.system("xplor -py xplor_build_link.py")
    #~ link = PDBModel("01_extended_link.pdb")
    
    link = Xtender(linkseq).run()

    ## Orientation
    link.xyz -= link.centerOfMass()
    link_end = link.xyz[N.where( link.maskFrom( 'name', ['OT1'] ))[0][0]] 
    link_start = link.xyz[N.where( link.maskFrom( 'name', ['HT2'] ))[0][0]] 
    orientVectors(link,link_end-link_start,direction)
    link_end = link.xyz[N.where( link.maskFrom( 'name', ['OT1'] ))[0][0]] 
    link_start = link.xyz[N.where( link.maskFrom( 'name', ['HT2'] ))[0][0]] 
    
    
    ## Linking... (1st step)
    
    a_last_C = a.xyz[N.where(a.maskFrom( 'name', ['C'] ))[0][-1]]
    b_first_N = b.xyz[N.where(b.maskFrom( 'name', ['N'] ))[0][0]]
    
    a.xyz += link_start - a_last_C
    b.xyz += link_end - b_first_N
    
    
    ## Fix joining sites planarity
 
    last_a_C = a.xyz[N.where(a.maskFrom( 'name', ['C'] ))[0][-1]]
    last_a_O = a.xyz[N.where(a.maskFrom( 'name', ['O'] ))[0][-1]]
    
    ini_link_N = link.xyz[N.where(link.maskFrom( 'name', ['N'] ))[0][0]]
    ini_link_Ca = link.xyz[N.where(link.maskFrom( 'name', ['CA'] ))[0][0]] 
    
    a,R = orientPlanes(modela = a, orig_plane = [last_a_O-last_a_C,ini_link_N-last_a_C], orig_perp = False, target_plane=[(last_a_C-ini_link_N),(ini_link_Ca-ini_link_N)],target_perp = False)
    
    
    ini_b_N = b.xyz[N.where(b.maskFrom( 'name', ['N'] ))[0][0]]
    ini_b_Ca = b.xyz[N.where(b.maskFrom( 'name', ['CA'] ))[0][0]] 
    
    last_link_C = link.xyz[N.where(link.maskFrom( 'name', ['C'] ))[0][-1]]
    last_link_O = link.xyz[N.where(link.maskFrom( 'name', ['O'] ))[0][-1]]
    
    b,R = orientPlanes(modela = b, orig_plane = [last_link_C-ini_b_N,ini_b_Ca-ini_b_N], orig_perp = False, target_plane=[(last_link_O-last_link_C) ,(ini_b_N-last_link_C)],target_perp = False)
       
    
    ## Linking together... (2nd step)

    a_last_C = a.xyz[N.where(a.maskFrom( 'name', ['C'] ))[0][-1]]
    b_first_N = b.xyz[N.where(b.maskFrom( 'name', ['N'] ))[0][0]]
    
    a.xyz += link_start - a_last_C
    b.xyz += link_end - b_first_N
    
    
    ## - Change chains and segments
    a.atoms['chain_id'] = ["A"]*len(a)
    b.atoms['chain_id'] = ["C"]*len(b)
    link.atoms['chain_id'] = ["B"]*len(link)
    
    link = cleanPdb(link)
    final = a.concat(link).concat(b)
    final.renumberResidues()
    final['serial_number'] = range( len(final) )
    final.atoms['segment_id'] = ["p1"]*len(final)
    final.report()
    final =  cleanPdb(final)
    
    fixed = []
    where = int(0.9*a.lenResidues())
    fixed += range(where,a.lenResidues())
    where = int(0.1*b.lenResidues())
    fixed += range(link.atoms['residue_number'][-1],where+link.atoms['residue_number'][-1]+1)
    fixed += range(link.atoms['residue_number'][0]-1,link.atoms['residue_number'][-1]+1)
    
        
    for i in range(len(final)):
        if final.atoms['residue_number'][i] in fixed:
            final.atoms['temperature_factor'][i] = 1.0
        else:
            final.atoms['temperature_factor'][i] = 0.0
    
    final.writePdb("initial.pdb")
    
    
    #~ os.system('vmd -dispdev text -e vmdscript')
    #~ os.system('namd2 namdscript')
    
    
    #~ final = PDBModel('in_min.pdb')
    #~ for i in range(len(final)):
        #~ if final.atoms['residue_number'][i] in fixed:
            #~ final.atoms['temperature_factor'][i] = 1.0
        #~ else:
            #~ final.atoms['temperature_factor'][i] = 0.0
    #~ final.writePdb("in_min.pdb") 
    
    #~ for i in range(0,3):
        #~ os.system('namd2 namdscript')
        #~ os.system('cp out_min.coor in_min.pdb')
    
    #~ os.system('cp out_min.coor final.pdb')


def joinProteins(a = None, b = None):
    """
    Binds covalently two proteins.
    
    @param a: First model to be linked.
    @type a: PDBModel
    @param b: Second model to be linked.
    @type b: PDBModel
    
    """
    
  
    
    cm1 = a.centerOfMass()
    
    cm2 = b.centerOfMass()
    
    a = cleanPdb(a)
    b = cleanPdb(b)
    
    
    ## Go to the origin!
    a.xyz -= cm1
    b.xyz -= cm2
    
    
    a.writePdb("initial.pdb")
    os.system('vmd -dispdev text -e vmdscript')
    a = PDBModel("in_min.pdb")
    a_end = a.xyz[N.where( a.maskFrom( 'name', ['OT1'] ))[0][0]] 
    a_start = a.xyz[N.where( a.maskFrom( 'name', ['HT2'] ))[0][0]] 
    a = cleanPdb(a)
    
    
    
    
    
    
    a_last_C = a.xyz[N.where(a.maskFrom( 'name', ['C'] ))[0][-1]]
    b_first_N = b.xyz[N.where(b.maskFrom( 'name', ['N'] ))[0][0]]
    
    b.xyz += a_end - b_first_N
    
    
    ## Fix joining sites planarity
 
    last_a_C = a.xyz[N.where(a.maskFrom( 'name', ['C'] ))[0][-1]]
    last_a_O = a.xyz[N.where(a.maskFrom( 'name', ['O'] ))[0][-1]]
    
    ini_b_N = b.xyz[N.where(b.maskFrom( 'name', ['N'] ))[0][0]]
    ini_b_Ca = b.xyz[N.where(b.maskFrom( 'name', ['CA'] ))[0][0]] 
    
    b,R = orientPlanes(modela = b, orig_plane = [last_a_C-ini_b_N,ini_b_Ca-ini_b_N], orig_perp = False, target_plane=[(last_a_O-last_a_C) ,(ini_b_N-last_a_C)],target_perp = False)
       
    
    ## Linking together... (2nd step)

    a_last_C = a.xyz[N.where(a.maskFrom( 'name', ['C'] ))[0][-1]]
    b_first_N = b.xyz[N.where(b.maskFrom( 'name', ['N'] ))[0][0]]
    
    b.xyz += a_last_C - b_first_N
    
    
    ## - Change chains and segments
    a.atoms['chain_id'] = ["A"]*len(a)
    b.atoms['chain_id'] = ["B"]*len(b)
   
    final = a.concat(b)
    final.renumberResidues()
    final['serial_number'] = range( len(final) )
    final.atoms['segment_id'] = ["p1"]*len(final)
    final.report()
    final =  cleanPdb(final)
    
    
    tata = a.lenResidues()
    tatb = b.lenResidues()
    
    free_residues = range(int(tata*0.8), int(tata+(tatb*0.2)))
        
  
    final.writePdb("initial.pdb")
    
    
    #~ os.system('vmd -dispdev text -e vmdscript')
    
    #~ final = PDBModel('in_min.pdb')
    #~ for i in range(len(final)):
        #~ if final.atoms['residue_number'][i] in free_residues:
            #~ final.atoms['temperature_factor'][i] = 0.0
        #~ else:
            #~ final.atoms['temperature_factor'][i] = 1.0
        
    #~ final.writePdb("in_min.pdb")
    
    forces = final.clone()
    for i in range(len(final)):
        forces.atoms['occupancy'][i] = 100.0
        if not forces.atoms['residue_number'][i] in free_residues:
            if forces.atoms['residue_number'][i]<tata:
                forces.xyz[i] = N.array([1000.,1000.,1000.])
            else:
                forces.xyz[i] = N.array([-1000.,-1000.,-1000.])
        else:
            forces.xyz[i] = N.array([0.,0.,0.])
    forces.writePdb('in_min_forces.pdb')
    
    os.system('namd2 namdscript')
    
    
    
    
    
    #~ for i in range(0,1):
        #~ os.system('namd2 namdscript')
        #~ os.system('cp out_min.coor in_min.pdb')
    
    os.system('cp out_min.coor final.pdb')
    return cleanPdb(PDBModel('final.pdb'))

if __name__ == '__main__':    
    #~ linkProteins(a = PDBModel('2AWT') , b= PDBModel('2CQJ'),distance = 70)
    #~ joinProteins(a = PDBModel('2AWT.pdb') , b= PDBModel('2CQJ.pdb'))
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
    
    p = cleanPdb(PDBModel('1JUN.pdb'))
    
    leu_A = p.takeChains([0])
    leu_B = p.takeChains([1])
    #~ linkProteins(a = PDBModel('1HUY.pdb') , b= leu_A, distance = 700)
    #~ os.system("cp initial.pdb A.pdb")
    #~ linkProteins(a = PDBModel('1HUY.pdb') , b= leu_B, distance = 700)
    #~ os.system("cp initial.pdb B.pdb")
    a = PDBModel('A.pdb')
    b = PDBModel('B.pdb')
    a_2 = a.takeChains([2])
    b_2 = b.takeChains([2])
    a.atoms['chain_index'] = ['A']*len(a)
    b.atoms['chain_index'] = ['B']*len(b)
    
    a_2.writePdb('lel.pdb')
    b_2.writePdb('lil.pdb')
    
    
    #recover original position of the leuzip
    
    
    av1 = leu_A.xyz[1]-leu_A.xyz[0]
    av2 = leu_A.xyz[2]-leu_A.xyz[0]
    
    a_2v1 = a_2.xyz[1]-a_2.xyz[0]
    a_2v2 = a_2.xyz[2]-a_2.xyz[0]
    
    
    
    newa,R = orientPlanes(modela = a, orig_plane = [a_2v1,a_2v2], orig_perp = False, target_plane=[av1,av2],target_perp = False)
    a_2 = newa.takeChains([2])
    
    bv1 = leu_B.xyz[1]-leu_B.xyz[0]
    bv2 = leu_B.xyz[2]-leu_B.xyz[0]
    
    b_2v1 = b_2.xyz[1]-b_2.xyz[0]
    b_2v2 = b_2.xyz[2]-b_2.xyz[0]
    
    newb,R = orientPlanes(modela = b, orig_plane = [b_2v1,b_2v2], orig_perp = False, target_plane=[bv1,bv2],target_perp = False)
    b_2 = newb.takeChains([2])
    
    cm_a = leu_A.centerOfMass()
    cm_b = leu_B.centerOfMass()
    
    newa.xyz += cm_a - a_2.centerOfMass()
    newb.xyz += cm_b - b_2.centerOfMass()
    
    newa.concat(newb).writePdb('result.pdb')
    newa.writePdb('result2.pdb')