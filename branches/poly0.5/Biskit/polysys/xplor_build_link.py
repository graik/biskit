import protocol, psfGen, pdbTool



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
    
seq = open('01_link_sequence.txt').read()
seq2 = single2longAA(seq)
psfGen.seqToPSF(seq2, segName='LINK', seqType='protein' )
protocol.genExtendedStructure()
#~ protocol.fixupCovalentGeom(maxIters=100,useVDW=1)
pdbTool.PDBTool("01_extended_link.pdb").write()