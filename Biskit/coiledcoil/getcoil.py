from coiledcoil import CoiledCoil
from choosecoil import CCStudy
import methods

CHAINS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def dataFileCreation( candidates_file = ""  ,target_seq = "",target_type = ("homodimer","parallel"), try_all = False):
    """
    Parses the candidates file and gets the master data file.
    
    @param candidates_file: Path for the candidates file.
    @type candidates_file: string
    
    @Candidates file format:
    The first line must be the path to a folder containing the pdb structures
    of our templates.
    
    Then it's followed by several lines with this syntax:
    
    PDB_FILE_NAME chainA [chainB]
    
    Where:
    
    PDB_FILE_NAME is a string with the file name of the pdb structure model we refer
    
    chainA is the chain to be used as reference for alignment
    
    and 
    
    chainB is the complementary chain (which is optional in the case we are
    using a homodimer)
    
    """
    try:
        file = open(T.dataRoot() + '/coiledcoil/'+table,"r")
        lineas = file.readlines()
    except IOError, msg:
        raise BiskitError('Cannot open score file %s.\n Error: %s' \
                          % (db, msg ) )
    file.close()
    
    assert (not(target_type[0]=="homodimer" or target_type[0]=="heterodimer")),"Options for target type are 'heterodimer' or 'homodimer'." 
    assert (target_seq == ""),"You have to define a target sequence."
    
    lineas = [ l.strip() for l in lineas ]
    basepath = lineas[0]
    
    
    if try_all:
        for 
    else:
        if target_type[1] == "parallel":
            ccdb = db or T.dataRoot() + '/coiledcoil/SOCKET_par_norm'
        elif:
            ccdb = db or T.dataRoot() + '/coiledcoil/SOCKET_antipar_norm'
        else:
            raise BiskitError("Options for target type are 'parallel' or 'antiparallel'." )
    
    file = open(basepath+"/_data","w")
    for l in lineas[1:]:
        line = l.split()
        new_line = line[0][:-3]+" struct:"+line[0]
        if 
        file.writeline(new_line)
    file.close()


def getCoilStructure ( type = ("parallel","homodimer"), db =""):
    """
    Function for coiled coil prediction.
    """
    assert (candidates_file != ""), " A candidates data file is needed for this function to work"
    
    if type[0] == "parallel":
        ccdb = db or T.dataRoot() + '/coiledcoil/SOCKET_par_norm'
    elif:
        ccdb = db or T.dataRoot() + '/coiledcoil/SOCKET_antipar_norm'
    else:
        ## ERROR
        pass
        
    
    study = CCStudy(cc)
    
    ## Data file creation (candidates parsing)




##############
## Test
##############
import Biskit.test as BT
import Biskit.tools as T
from alignment import PirAlignment

class Test(BT.BiskitTest):
    """ Test cases for Coiled Coil Utils"""

    def prepare(self):
        pass

    def cleanUp( self ):
        pass
        
    def test_getCoil(self):
        """doStudy function test case"""
        
        
if __name__ == '__main__':
    BT.localTest()    
       