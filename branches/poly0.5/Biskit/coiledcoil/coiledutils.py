from string import letters,digits
from numpy import array
import Biskit.molUtils as MU

def scoreAlignFun(a,b):
    printable = letters + digits
    return  len(printable) - abs(printable.find(a)-printable.find(b))
    
    
def scores2String(a=[],b=[],reversed = False):
    printable = letters + digits
    ar = "";br = ""
    maxp = max(a+b)
    na = array(a);nb = array(b)
    if reversed:
        na = maxp-na;nb = maxp - nb
    maxp = max(a+b)
    na = na/maxp;nb = nb/maxp
    maxi = len(printable)-1

    for i in na:
        ar += printable[int(maxi*(i))]

    for i in nb:
        br += printable[int(maxi*(i))]

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

def sameRegister(hept1 = "",hept2 ="",chain=""):
    """
    This function returns True when two heptads (hept1,hept2) have
    the same registers in chain (so their registers have to be 'abcdefg'
    for both.
    
    @param hept1: First heptad to compare.
    @type hept1: string
    @param hept2: Second heptad to compare.
    @type hept2: string
    @param chain: String to which this heptads belong.
    @type chain: string
    
    @return: True if both share registers.
    @rtype: boolean
    """
    a = chain.find(hept1)
    b = chain.find(hept2)
    
    return (a-b)%7 == 0 


def alignd(a="",b="",dict={}):
    
    acc = []
    for pos in range(len(a)-len(b)+1):
        aux = 0
        for i in range(len(b)):
            aux += float(dict[(a[pos+i],b[i])])
        acc.append((aux,pos))
    return sorted(acc,reverse = True)

def alignf(a="",b="",fun=None):
    
    acc = []
    for pos in range(len(a)-len(b)+1):
        aux = 0
        for i in range(len(b)):
            aux += float(fun(a[pos+i],b[i]))
        acc.append((aux,pos))
    return sorted(acc,reverse = True)


def flatten ( chain = "", register = "", mycc = None,window_length = 7):
    """
    This function returns a list with the scores of each heptad in chain.
    Incomplete heptads are not used.
    
    @param chain: Chain to be flattened.
    @type chain: string
    @param register: Register ofthe chain.
    @type heptad: string
    
    @return: A list with the scores of each heptad in chain.
    @rtype: list (float)
    """
    
    start = register.find("a")
    end  = start + len(chain[start:]) - len(chain[start:])%window_length
  
    
    scores = []
    pos = 0
    ac = 0 

    for i in chain[start:end]:
        ac+= mycc.scores[MU.single2longAA(i)[0]][pos] or 0
        pos = (pos+1)%window_length
        if pos == 0:
            scores.append(ac)
            ac = 0
    
    return scores


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
    def test_alignment(self):
        """ testing of mini-alignment functions"""
        from coiledalign import CoiledAlign
        ca = CoiledAlign()
        
        res = alignd("RRRLLLLLLLRRRRRRLLLLLLLRRR", "LLLRRRR",ca.like_scores)
        self.assertEqual(res[0],(32.0, 7))
        
        res =  alignd("RRRRRRL","LLLRRRR",ca.charge_scores)
        self.assertEqual(res[0],(2.0, 0))
        
        res = alignf("RRRLLLLLLLRRRRRRLLLLLLLRRR", "LLLRRRR",scoreAlignFun)
        self.assertEqual(res[0],(434.0, 7))
        
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
            print
            print scores2String([0.1,0.3,0.5,0.6],[1,1.5])
        self.assertEqual(scores2String([0.1,0.3,0.5,0.6],[1,1.5]),('emuy', 'O9'))
        self.assertEqual(scores2String([0.1,0.3,0.5,0.6],[1,1.5],True),('4WOK', 'ua'))

if __name__ == '__main__':
    BT.localTest()    
    