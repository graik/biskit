

def lendepth(this):
    len = 0
    
    try:
        for i in this.keys():
            len += lendepth(this[i]) 
        return len
    except:
        try:
            for i in this:
                len += lendepth(i) 
            return len
        except:
            return 1
    
    
##############
## Test
##############
import Biskit.test as BT
import Biskit.tools as T


class Test(BT.BiskitTest):
    """ Test cases for fold data creation"""

    def prepare(self):
        pass

    def cleanUp( self ):
        pass
    
    def test_lendepth(self):
        """Lendepth test case"""
        a1 = range(4)
        a2 = range(5)
        a3 = range(5)
        
        b1 = range(3)
        b2 = {"ab":1,"bb":2}

        a = [a1,a2,a3]
        b = [b1,b2]
        
        d = {"a":a,"b":b}
        
        self.assertEqual( lendepth(d),19)
        
if __name__ == '__main__':
    BT.localTest()    
 