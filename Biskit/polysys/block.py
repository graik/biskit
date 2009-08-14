from interval import Interval


class BlockEntity:
    """
    Wrapper for diferent kind of classes that can be used for data storage. It has some common functions
    like the ones used to store data.
    Different blocks can have the same main entity.
    """
    
    def __init__(self,name= "Label",structure = None):
        """
        Creates a BlockEntity.
        
        @param name: Is the identifier of a BlockEntity. Must be unique.
        @type name: string
        """
        self.name = name
        self.structure = structure
        
    
    def __str__(self):
        """
        Returns a string representation of a BlockEntity. Its description is basically "Entity" followed by its name.
        
        @return : String representation.
        @rtype: string
        """
        return "[Entity: "+self.name+"]"
        
    def onInsertion(self,myblock = None,myassembly= None):
        pass
    
    def run(self):
        pass
        
    def getStructData(self,atoms):
        pass

class SimpleBlockEntity (BlockEntity):
    def __init__(self,name= "Label",structure = None):
        """
        Creates a BlockEntity.
        
        @param name: Is the identifier of a BlockEntity. Must be unique.
        @type name: string
        """
        
        BlockEntity.__init__(self,name,structure)
    
    
    def run(self):
        return self.structure
    
    def getStructData(self,atoms):    
        return self.take(atoms)
    
class Block(object):
    """
    Each Block is a tree-like structure identified by its name (which must be unique).
    If it has a father block, it has to share its entity wih it.
    """
    
    def __init__(self,name, father = None,entity = None ):
        """
        Creates a Block.
        @param name: Is the identifier of a Block. Must be unique.
        @type name: string
        @param entity: Entity which this block is going to contain. Can be None.
        @type entity: BlockEntity
        @param father: Father Block in a tree-like structure. Can be None.
        @type father: Block
        """
        self.name = name
        
        
        
        if father:
            assert isinstance(father,Block) ,'father is not a Block'
        
        self.entity = entity
        self.intervals = []
        self.sons = []
        
        self.father = father
        if self.father != None :
            self.father.sons.append(self)
        
        if self.father and not self.entity:
            self.entity = self.father.entity
    
    def __dettachFromFather(self):
        """
        Removes son from father list and sets it to None.
        """
        if self.father:
            self.father.sons.remove ( self )
            self.father = None
    
    def setFather (self,f):
        """
        Changes the father of the current block in a consistent way.
        
        @param f: The new father.
        @type f: Block
        """
        
        if self.father:
            self.__dettachFromFather()
            
        self.father = f
        self.father.sons.append(self)
        self.entity = self.father.entity
        
    def addInterval(self, s,e=-1):
        """
        Function for adding an interval from the interval list.
        Can be called with only one parameter, then 's' is an interval instance.
        
        @param s: Interval to be removed or starting point of the interval.
        @type s: Interval or int
        @param e: Ending point of the interval.
        @type e: int
        """
        
        if e==-1:
            assert not isinstance(self,Interval) ,'Please check parameter types.'
            i = s
        else:
            i = Interval(s,e)
        
        assert( not i.checkOverlap(self.intervals)), "Interval overlap for the same Block"
    
        self.intervals.append(i)	
        
    def remInterval(self, s,e = -1):
        """
        Function for removing an interval from the interval list.
        Can be called with only one parameter, then 's' is an interval instance.
        
        @param s: Interval to be removed or starting point of the interval.
        @type s: Interval or int
        @param e: Ending point of the interval.
        @type e: int
        """
        if e==-1:
            assert not isinstance(self,Interval) ,'Please check parameter types.'
            i1 = s
        else:
            i1 = Interval(s,e)
            
        for i in self.intervals:
            if i == i1:
                self.intervals.remove(i1)
                return
        
    def __cmp__(self,o):
        """
        Compares two blocks. Two blocks are the same if they have the same name. 
        """
         
        if o != None and self.name == o.name :
            return 0
        else:
            return -1
    
    def __str__(self):
        """
        Returns a string representation of a Block. Consists of the name of the block, its intervals, and
        entity description.
        
        @return : String representation.
        @rtype: string
        """
        mystring = "[Block "+self.name+"]"
        
        if self.intervals != []:
                for i in self.intervals:
                    mystring += "\n"+ i.__str__()
        
        if self.entity != None :
            mystring += "\n"+self.entity.__str__()
            
        return mystring
        
    def getStringTree(self,level = 0):
        """
        Creates a string containing a treelike representation of this block and subblocks.
        
        @param level: Used to keep track of the indentation. Calling the function with other values 
                will only affect indentation.
        @type level: int
        """
        mystring = self.name
        
        nle = level+1
        
        for i in self.sons :
            mystring += '\n'
            mystring += "--"* level
            mystring +="-->SubBlock : "
            mystring +=i.getStringTree(nle) 
        
        return mystring
        
    def getModel (self) :
        """
        Function for retrieving structural data from entities.
        """
        r  = []
        if self.entity:
            for i in self.intervals:
                r += i.getRange()
            
            self.entity.getStructData(r)

    def onInsertion(self,assembly= None):
        if self.father == None: # so is only executed for the "top" block
            if self.entity:
                self.entity.onInsertion(myassembly = assembly)
            
##############
## Test
##############
import Biskit.test as BT
import Biskit.tools as T
from Biskit.PDBModel import PDBModel

class Test(BT.BiskitTest):
    """ Test cases for Polyfret"""
    
    def prepare(self):
        pass

    def cleanUp( self ):
        pass
    
    def test_Blocks(self):
        """Block test cases"""
        
        b1 = Block("my_level1_block")
        b11 = Block("my_level2_block_a",father = b1)
        b12 = Block("my_level2_block_b",father = b1)
        b111 = Block("my_level3_block_a",father = b11)
        b1111 = Block("my_level4_block_a",father = b111)
        b112 = Block("my_level3_block_b",father = b11)
        b121 = Block("my_level3_block_a",father = b12)
        
        if self.local:
            print b1
            print b1.getStringTree()
        
        self.assertEqual(b121.father.father.name,"my_level1_block")
        
        b111.setFather(b1)
        
        self.assertEqual(b111.father.name,"my_level1_block")
        self.assertEqual(b111 in b111.father.sons,True)
        
        if self.local:
            print b1.getStringTree()
        
        b1.addInterval(Interval(1,4))
        b1.addInterval(Interval(7,10))
        self.assertEqual(len(b1.intervals),2)
        
        b1.getModel()
        
if __name__ == '__main__':
    BT.localTest()