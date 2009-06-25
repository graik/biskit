from coiledcoil import CoiledCoil
from choosecoil import CCStudy

CHAINS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def dataFileCreation( candidates_file = ""  ):
    """
    Parses the candidates file and gets the master data file.
    
    @param candidates_file:
    @type candidates_file: 
    
    Candidates file format:
    
    """
    try:
        lineas = open(T.dataRoot() + '/coiledcoil/'+table,"r").readlines()
    except IOError, msg:
        raise BiskitError('cannot open score file %s.\n Error: %s' \
                          % (db, msg ) )
        


"""
Function for coiled coil prediction.
"""

def getCoilStructure ( , type = ("parallel","homodimer"), db =""):
    assert (candidates_file != ""), " A candidates data file is needed for this function to work"
    
    if type[0] == "parallel":
        cc = db or T.dataRoot() + '/coiledcoil/SOCKET_par_norm'
    elif:
        cc = db or T.dataRoot() + '/coiledcoil/SOCKET_antipar_norm'
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
       