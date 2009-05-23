from string import letters,digits
from numpy import array

def scoreAlignFun(a,b):
    return abs(ord(a)-ord(b) / max(ord(a),ord(b)))
    
    
def scores2String(a,b):
    printable = letters + digits
    
    ar = "";br = ""
    maxp = max(max(a),max(b))
    na = array(a);nb = array(b)
    na = na / maxp;nb = nb / maxp
    maxi = len(printable)-1
    for i in na:
        ar += printable[int(maxi*i)]
    for i in nb:
        br += printable[int(maxi*i)]
    return ar,br
    
    
def getRegister(heptad="",chain="",window_length = 7):
        
    index = chain.find(heptad)
    if index % window_length == 0:
        head =""
    else:
        head = "abcdefg"[-(index % window_length):]
    reg1 = head+"abcdefg"*((len(chain)-(index % window_length))/7) 
    reg2 = reg1+"abcdefg"[:len(chain)-len(reg1) ]
    
    return reg2    


##############
## Test
##############
import Biskit.test as BT
import Biskit.tools as T

class Test(BT.BiskitTest):
    """ Test cases for Coiled Coil Utils"""

    def prepare(self):
        pass

    def cleanUp( self ):
        pass
        
    def test_getRegister(self):
        """getRegister function test case"""
        chain="LEIRAAFLRRRNTALRTRVAELRQRVQRLRNIVSQYETRYGPL"
        reg = getRegister("LRTRVAE",chain)
        
        if self.local:
            print
            print chain
            print reg
        
        self.assertEqual(reg,"abcdefgabcdefgabcdefgabcdefgabcdefgabcdefga")
        
        reg = getRegister("VSQYETR",chain)
        
        if self.local:
            print chain
            print reg
        self.assertEqual(reg,"defgabcdefgabcdefgabcdefgabcdefgabcdefgabcd")
    
    def test_scores2String(self):
        """scores2String function test case"""
        if self.local:
            print scores2String([0.1,0.3,0.5,0.6],[1,1.5])
        self.assertEqual(scores2String([0.1,0.3,0.5,0.6],[1,1.5]),('emuy', 'O9'))

if __name__ == '__main__':
    BT.localTest()    
    