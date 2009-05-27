


class Alignment :
    
    def __init__ (self, chain_a = "",chain_b="",pos = 0):
        self.chain_a = chain_a
        self.chain_b = chain_b
        self.pos = pos
        self.__padding()
    
    def __padding (self):
        len_a =  len( self.chain_a )
        len_b =  len( self.chain_b )
        
        if self.pos < 0:
            self.chain_a = "-"*(abs(self.pos))+self.chain_a
            self.pos = 0
            len_a =  len( self.chain_a )
        
        if len_b + self.pos > len_a:
            self.chain_a += "-"*(len_b + self.pos - len_a)
            len_a =  len( self.chain_a )
            
        if len_a - (len_b + self.pos) > 0:
            self.chain_b += "-"*(len_a - (len_b + self.pos))
        
        self.chain_b = ("-"*(self.pos))+self.chain_b
    
    def crop (self,what = "",where = 50):
        chains = []
        len_c = len(what)
        
        i = 0
        for i in range(len_c / where):
            chains.append(what[(i*where):((i+1)*where) ])

        if len_c / where == 0:
            chains.append(what)
        else:
            chains.append(what[((i+1)*where):])

        return chains
        
    def writePir (self,path = "",prot_name='XXXX'):
        file = open(path,"w")
        
        file.write(">P1;"+prot_name+"\n")
        file.write("structure:"+prot_name+": : : : : : : :\n")
        file.writelines(self.crop(self.chain_a))
        file.write("\n*\n>P1;target\n")
        file.write("sequence:target: : : : : : : :\n")
        file.writelines(self.crop(self.chain_b))
        file.write("\n*")
        file.close()
        
##############
## Test
##############
import Biskit.test as BT
import Biskit.tools as T

class Test(BT.BiskitTest):
    """ Test cases for Alignment"""

    def prepare(self):
        pass

    def cleanUp( self ):
        pass
        
    def test_align_and_pad(self):
        """Align test cases"""
        a = Alignment ( "1234567890","456",0)
        self.assertEqual(a.chain_b,"456-------")
        if self.local:
            print 
            print a.chain_a
            print a.chain_b
        a = Alignment ( "1234567890","456",3)
        self.assertEqual(a.chain_b,"---456----")
        if self.local:
            print
            print a.chain_a
            print a.chain_b
        a = Alignment ( "1234567890","456",8)
        self.assertEqual(a.chain_b,"--------456")
        if self.local:
            print
            print a.chain_a
            print a.chain_b
        a = Alignment ( "1234567890","456",-3)
        self.assertEqual(a.chain_b,"456----------")
        if self.local:
            print
            print a.chain_a
            print a.chain_b
        
    def test_crop(self):
        """Crop test"""
        a = Alignment ( )
        self.assertEqual(len(a.crop("1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456")),3)
        self.assertEqual(len(a.crop("1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456")[2]),6)
        self.assertEqual(len(a.crop("123456789012345678901234")),1)
        
        
    def test_write(self):
        """Write test"""
        a = Alignment ( "1234567890","456",0)
        a.writePir("lol","lalala")
        
if __name__ == '__main__':
    BT.localTest() 