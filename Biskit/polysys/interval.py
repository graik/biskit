from polysyserrorhandler import PolySysErrorHandler


class Interval:
    """
    Class for interval storage. An interval is the definition of a range.
    """
    def __init__(self,s,e,verbose = True):
        """
        Interval creation. Interval goes from 's' to 'e'. An interval must be defined by two 2 @positive numbers
        and the first one shall be smaller than the second (but if this happen it will automatically swap the 
        ending and starting points).
        
        @param s: Starting point of an interval (range).
        @type s: int
        @param e: Ending point of an interval.
        @type e: int
        """
        self.start = s
        self.end  = e
        self.ehandler = PolySysErrorHandler(verbose)
        self._Interval__checkConsistence()
        
    
    def __checkConsistence(self):
        """
        Checks if the interval has been created following some format rules.
        """
        if( self.start<0 or self.end<0 ) :
            self.ehandler.error( "Start or end can't be negative numbers when defining the range.")
            self.start = 0
            self.end  = 0
        if(self.start>self.end):
            i = self.start
            self.start=self.end
            self.end=i
            self.ehandler.warning("["+self.__str__()+"] Swapping start and end.")
            
    def checkOverlap(self, others = []):
        """
        Checks if this interval overlaps with the others in param "others".
        
        param others: All the intervals to check overlapping.
        type others: Interval list
        
        return: Will return 'True' if any Interval in 'others' overlaps with the caller.
        rtype: bool
        """
        for i in others:
            if not ( (self.end > i.end and self.start > i.end) or (self.start<i.start and self.end < i.start) ):
                return True
        
        return False
        
    def __len__(self):
        return self.end +1 - self.start
    
    def __str__(self):
        """
        Returns a string representation of an Interval. Consists of the range if the interval.
        
        @return : String representation.
        @rtype: string
        """
        return '[Start: %d , End: %d]'%(self.start, self.end)
        
    def __cmp__(self,o):
        """
        Interval comparison.
        
        @return: True if the two intervals have the same start and ending points.
        @rtype: bool
        """
        if isinstance(o,Interval):
            if o.start == self.start and o.end == self.end :
                return 0
            else :
                return -1
        elif isinstance(o,Spots):
            return cmp(o,self)
        else:
            return -1
    
    def getRange( self ):
        """
        Creates a list containing the vlaues of the range defined by start and end (where end is included!)
        """
        return range(self.start,self.end+1)
        
    def shift(x):
        """
        Shiftes the interval an x ammount.
        """
        self.start = self.start + x
        self.end = self.end + x
        

class Spots :
    """
    Spots is a class for describing properties wich can be defined in lots of sectors of the structure. It can also
    be used as an interval, or even as a group of intervals, but it's not recommended.
    """
    
    def __init__(self, mylist = []):
        """
        Creation function.
        
        @param mylist: list containing an ensemble of spots for the propertie (for instance the index of all
                    Ca carbons.
        @type mylst: int list
        """
        self.desclist = mylist
        
    def __cmp__(self,o):
        """
        Spot comparison.
        
        @return: True if the two spots are exactly the same.
        @rtype: bool
        """
        if isinstance(o,Spots):
            if self.desclist == o.desclist :
                return 0
            else :
                return -1
        elif isinstance(o,Interval):
            if len(self) == len(o) and self.desclist[0] == o.start:
                return 0
            else:
                return -1
        else:
            return -1
            
    def __len__(self):
        return len(self.desclist)
    
    def __str__(self):
        """
        Returns a string representation of a Spots.
        
        @return : String representation.
        @rtype: string
        """
        return str(self.desclist)+"\n"
        
    def getRange( self ):
        """
        Returns the list
        """
        return self.desclist

    def shift(x):
        """
        Shiftes the interval an x ammount.
        """
        for i in range(len(self.desclist)):
            self.desclist[i] = self.desclist[i] + x

##############
## Test
##############
import Biskit.test as BT
import Biskit.tools as T

class Test(BT.BiskitTest):
    """ Test cases for Interval"""
    
    def prepare(self):
        pass

    def cleanUp( self ):
        pass
    
    def test_Interval(self):
        """Interval test cases """
        
        # Overlapping
        a = Interval(1,3)
        b = Interval(6,8)
        c = Interval(4,5)
        d1 = Interval(2,5)
        e = Interval(5,9)
        m = Interval(1,4)
        n = Interval(9,10)
        self.assertEqual (  c.checkOverlap([m,n]), True )
        self.assertEqual (  c.checkOverlap([a,b]), False )
        self.assertEqual (  c.checkOverlap([b]), False )
        self.assertEqual (  c.checkOverlap([d1]),True )
        self.assertEqual (  c.checkOverlap([a,b,d1]),True )
        self.assertEqual (  e.checkOverlap([b]),True )
        self.assertEqual (  b.checkOverlap([e]),True )
        self.assertEqual (  e.checkOverlap([a]),False )
        
        # Consistence 
        try:
            a = Interval(-1,3)
        except:
            pass
        else:
            self.assertEqual (1,2)
            
        a = Interval(3,1,False)
        self.assertEqual ( "Swapping start and end" in a.ehandler.lastWarning,True)
        
        
    def test_Spots(self):
        """Spots test cases """
        
        s1 = Spots(range(0,31))
        s2 = Spots(range(0,31))
        s3 = Spots(range(0,6))
        c = Interval(0,5)
        self.assertEqual (   len(s1), 31)
        self.assertEqual (  s1, s2)
        self.assertNotEqual (  s1, s3)
        self.assertNotEqual (  s1 ,c)
        self.assertNotEqual (  c , s1)
        self.assertNotEqual (  c , s2)
        self.assertEqual (  c , s3)
        self.assertEqual (  s3, c)
        
if __name__ == '__main__':

    BT.localTest()