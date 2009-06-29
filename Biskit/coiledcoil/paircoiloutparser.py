class PairCoilResult:
    """
    Class for socket results storage.
    """
    def __init__(self):
        """
        Instantiation.
        A socket result just has a sequence and its register.
        """
        self.chain = ""
        self.register = ""
        
    def __str__ (self):
        """
        Sring conversion and printing handling.
        """
        mystr = "--- Result ---\n"
        mystr += "Chain:    "+str(self.chain)+"\n"
        mystr += "Register: "+str(self.register)
        return mystr
        

def parse(path):
    """
    Parsing for PairCoil2 regular output files.
    
    @return: Result of the parsing.
    @type: list {PairCoilResult}
    """
    
    file = open(path,"r")
    
    lines = file.readlines()
    
    file.close()
    
    
    lineas = [ l.strip() for l in lines ]
    
    
    results = {}
    processing_seq = False
    processing_reg = False
    
    result = PairCoilResult()
    
    for l in lineas:
        contents = l.split()
        
        if processing_seq:
            contents = l.split()
            if len(contents) > 1 :
                result.chain +=contents[1]
        
        if processing_reg:
            contents = l.split()
            if len(contents) > 1 and contents[0][0] != '[':
                result.register +=contents[1]
        
        if '0123456789' in l:
            if processing_seq == False:
                processing_seq = True
            else:
                processing_seq = False
                processing_reg = True
                
    return result
##############
## Test
##############
import Biskit.test as BT
import Biskit.tools as T

class Test(BT.BiskitTest):
    """ Test cases for CoiledCoil"""

    def prepare(self):
        pass

    def cleanUp( self ):
        pass
        
    def test_general(self):
        """General parsing"""
        r1 =  parse('/home/victor/poly0.5/Biskit/testdata/coiledcoil/paircoutex1')
        r2 =  parse('/home/victor/poly0.5/Biskit/testdata/coiledcoil/paircoutex2')
        
        if self.local:    
            print r1
            print r2
        
        self.assertEqual(r2.chain,'GSHMPLLSIARQEEEMKEQLKQMDKMKEDLAKTERIKKELEEQNVTLLEQKNDLFGSMKQLEDKVEELLSKNYHLENEVARLKKLVGER')
        self.assertEqual(r2.register,'.....efgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcd')
        self.assertEqual(r1.chain,'RMKQLEDKVEELLSKKYHLENEVARLKKLVGER')
        self.assertEqual(r1.register,'gabcdefgabcdefgabcdefgabcdefgabcd')
        
        

        
if __name__ == '__main__':
    BT.localTest()    


        