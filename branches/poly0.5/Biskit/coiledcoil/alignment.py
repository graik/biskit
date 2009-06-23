


class PirAlignment :
    """ 
    Class for generating pairwise multichain alignments
    """
    
    def __init__ (self, chains=[("","")],pos = [0]):
        self.chains = chains
        self.processed_a = ""
        self.processed_b = ""
        self.fasta_b = ""
        for c in range(len(self.chains)):
            self.chain_a = self.chains[c][0]
            self.chain_b = self.chains[c][1]
            self.pos = pos[c]
            self.__padding()
            self.processed_a += self.chain_a
            self.processed_b += self.chain_b
            self.fasta_b = += self.chain_b
            if len(self.chains)-1!= c:
                self.processed_a += "/"
                self.processed_b += "/"
    
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
            chains.append(what[(i*where):((i+1)*where) ]+"\n")

        if len_c / where == 0:
            chains.append(what+"\n")
        else:
            chains.append(what[((i+1)*where):]+"\n")

        return chains
        
    def writePir (self,path = "",prot_name='XXXX'):
        file = open(path,"w")
        
        file.write(">P1;"+prot_name+"\n")
        file.write("structure:"+prot_name+": : : : : : : :\n")
        file.writelines(self.crop(self.processed_a))
        file.write("*\n>P1;target\n")
        file.write("sequence:target: : : : : : : :\n")
        file.writelines(self.crop(self.processed_b))
        file.write("*")
        file.close()
        
    def writeFasta (self,path = "",prot_name='XXXX'):
        file = open(path,"w")
        
        file.write(">"+prot_name+"\n")
        
        file.writelines(self.crop(self.fasta_b,60))
        
        file.close()
        
    def __str__(self):
        prot_name = "XXXX"
        str = ""
        str += ">P1;"+prot_name+"\n"
        str += "structure:"+prot_name+": : : : : : : :\n"
        for s in self.crop(self.processed_a):
            str += s
        str += "*\n>P1;target\n"
        str += "sequence:target: : : : : : : :\n"
        for s in self.crop(self.processed_b):
            str += s
        str += "*"
        return str
        
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
        a = PirAlignment ( [("1234567890","456")],[0])
        self.assertEqual(a.chain_b,"456-------")
        if self.local:
            print 
            print a.processed_a
            print a.processed_b
        a = PirAlignment ( [("1234567890","456")],[3])
        self.assertEqual(a.chain_b,"---456----")
        if self.local:
            print
            print a.processed_a
            print a.processed_b
        a = PirAlignment ( [("1234567890","456")],[8])
        self.assertEqual(a.chain_b,"--------456")
        if self.local:
            print
            print a.processed_a
            print a.processed_b
        a = PirAlignment ( [("1234567890","456")],[-3])
        self.assertEqual(a.chain_b,"456----------")
        if self.local:
            print
            print a.processed_a
            print a.processed_b
        
    def test_crop(self):
        """Crop test"""
        a = PirAlignment ( )
        self.assertEqual(len(a.crop("1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456")),3)
        self.assertEqual(len(a.crop("1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456")[2]),7)
        self.assertEqual(len(a.crop("123456789012345678901234")),1)
        
        
    def test_write(self):
        """Write test"""
        a = PirAlignment ( [("1234567890","456")],[0])
        a.writePir(T.testRoot()+"/coiledcoil/short.pir","MYPR")
        ## Long multi-alignment
        chain1a ="1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456"
        chain2a = "abcdefghijk"
        chain1b = "3456789"
        chain2b = "ghij"
        pos = [23,6]
        a = PirAlignment([(chain1a,chain2a),(chain1b,chain2b)],pos)
        a.writePir(T.testRoot()+"/coiledcoil/long.pir","MYPR")
        
if __name__ == '__main__':
    BT.localTest() 