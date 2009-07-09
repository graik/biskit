def pad(n,k=7,padchar= " "):
    
    if isinstance(n,int):
        number = "%-d"%(n)
    else:
        number = "%-.3f"%(n)
    number = padchar*(k-len(number))+number
    return number     

def cutdepth(this,thres):
    
    try:
       
        for i in this.keys():
            this[i] = cutdepth(this[i],thres) 
    except:
        try:
            
            for i in range(len(this)):
                this[i] = cutdepth(this[i],thres) 
        except:
            
            if abs(this) < thres:
                return 0
            else:
                return this
    return this    
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
    
    def test_pad(self):
        """Pad test case"""   
        self.assertEqual( pad(4.5),"  4.500")
    
    def test_cutdepth(self):
        """Threshold cut test case"""
        a = [[0.00002,0.1,-0.0002],[1e-23,3]]
        print 
        print a
        print cutdepth(a,0.001)
if __name__ == '__main__':
    BT.localTest()    
 