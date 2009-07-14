import psfGen, protocol, pdbTool
import random as rand
import protocol

aaDicStandard =\
              {'asp':'D', 'glu':'E', 'lys':'K', 'his':'H', 'arg':'R',
               'gln':'Q', 'asn':'N', 'ser':'S', 'asx':'B', 'glx':'Z',
               'phe':'F', 'trp':'W', 'tyr':'Y',
               'gly':'G', 'ala':'A', 'ile':'I', 'leu':'L', 'cys':'C',
               'met':'M', 'thr':'T', 'val':'V', 'pro':'P' }

def single2longAA( seq ):
    """
    Convert string of 1-letter AA code into list of 3-letter AA codes.
    
    @param seq: amino acid sequence in 1-letter code
    @type  seq: str
    
    @return: list with the amino acids in 3-letter code
    @rtype: [str]
    """
    ## invert AA dict
    invTab = {}

    for key in aaDicStandard:
        invTab[ aaDicStandard[key] ] = key

    result = []
    for aa in seq:
       
        aa = aa.upper()
        result += [ invTab[aa].upper() ]
       
    return result


protocol.initTopology(('protein'))
protocol.initParams(('protein'))

seq = open('01_sequence.txt').read()
seq2 = single2longAA(seq)
resid_a = open('01_resid_a.txt').read()
resid_b = open('01_resid_b.txt').read()
resid_link = open('01_resid_link.txt').read()


psfGen.seqToPSF(seq2,seqType='prot',startResid=1, segName='p1')

xplor.command( 'vector do (x = 0) (all)' )
xplor.command( 'vector do (y = %4.2f) (all)' % rand.random() )
xplor.command( 'vector do (z = %4.2f) (all)' % rand.random() )

## protocol.genExtendedStructure(maxFixupIters=50)

## read structures
protocol.initCoords( ['01_domain1.pdb', '01_domain2.pdb'] )

print '(attr x=0.0) and %s'%resid_a
protocol.fixupCovalentGeom(sel=AtomSel('(attr x=0.0) and %s'%resid_a),
                           maxIters=20, verbose=1 )
protocol.fixupCovalentGeom(sel=AtomSel('(attr x=0.0) and %s'%resid_b),
                           maxIters=20, verbose=1 )

## fix linker
xplor.command( 'vector do (x = %4.2f) (%s)' % (rand.random(),resid_link))
xplor.command( 'vector do (y = %4.2f) (%s)' % (rand.random(),resid_link))
xplor.command( 'vector do (z = %4.2f) (%)' % (rand.random(),resid_link))

protocol.fixupCovalentGeom(sel=AtomSel('%s'%resid_link), maxIters=20)

pdbTool.PDBTool("aux_out.pdb").write()
