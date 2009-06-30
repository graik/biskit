class SocketResult:
    """
    Class for socket results storage.
    For output file extraction. It stores:
    - Socket Id for the coiled coil.
    - Sense of the coiled coil (parallel or antiparallel)
    - Oligomerization type (2 for dimers, 3 for trimers...)
    - Chain sequences for the coiled coils
    - Registers predicted for these sequences.
    """
    
    def __init__(self):
        """
        Instantiation function.
        """
        self.coilId = 0
        self.sense = "parallel"
        self.chains = {}
        self.registers = {} 
        self.oligomerization = 2
        
    def __str__ (self):
        """
        Sring conversion and printing handling.
        """
        mystr = "--- Result ---\n"
        mystr += "ID: " + str(self.coilId)+"\n"
        mystr += "Oligo: " + str(self.oligomerization)+"\n"
        mystr += "Sense: " +str(self.sense)+"\n"
        mystr += "Chains: "+str(self.chains)+"\n"
        mystr += "Register: "+str(self.registers)
        return mystr

def parse (path):
    """
    Regular "non detailed" socket output files parsing.
    It just runs along each line detecting info chunks and then extracting
    info to "results".
    
    @return: List of results of coiled coils predicted by socket (see @SocketResult).
    @type: list {SocketResult}
    """
    
    file = open(path,"r")
    
    lines = file.readlines()
    
    file.close()
    
    
    lineas = [ l.strip() for l in lines ]
    
    
    results = {}
    coils = {}
    for l in lineas:
        contents = l.split()
        
        if len(contents)>4 and contents[0] == "coiled" and contents[1] == "coil" and not 'subset' in l:
            
            result = SocketResult()
            
            result.coilId = contents[2][:-1]
            
            coils[contents[2][:-1]] = []
            
            result.oligomerization =  int(contents[3])
            
            where = 0
            for i in range(len(contents)):
                if contents[i] == "frequency":
                    where = i 
            
            for i in range(5,where):
                result.chains[contents[i]] = ""
                result.registers[contents[i]] = ""
                coils[contents[2][:-1]].append(contents[i])
                
            results[result.coilId] = result
        
        
        if len(contents)>3 and contents[2] == "coiled" and contents[3] == "coil" :
            if not contents[4] == "packing:":
                results[contents[5]].sense = contents[6][1:]
                
        if len(contents)>2 and contents[0] == "assigning" and contents[1] == "heptad":
            helix = contents[4]
            for c in coils.keys():
                if helix in coils[c]:
                    mycoil = c
            
        if len(contents)>1 and contents[0] == "sequence":
            mysequence = contents[1]
            
        if len(contents)>1 and contents[0] == "register":
            myregister = l[9:]
            ## Intersection
            inter = ""
            for k in range(len(myregister)):
                if myregister[k] != " ":
                    inter += mysequence[k]
            results[mycoil].chains[helix] = inter
            results[mycoil].registers[helix] = myregister.split()[0]
            
    return results


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
        r1 =  parse('/home/victor/poly0.5/Biskit/testdata/coiledcoil/socketoutex1')
        r2 =  parse('/home/victor/poly0.5/Biskit/testdata/coiledcoil/socketoutex2')
        r3 =  parse('/home/victor/poly0.5/Biskit/testdata/coiledcoil/socketoutex3')
        r4 =  parse('/home/victor/poly0.5/Biskit/testdata/coiledcoil/socketoutex4')
        
        
        if self.local:
            for o in r1:
                print r1[o] 
            print
            
            for o in r2:
                print r2[o] 
            print
            
            for o in r3:
                print r3[o] 
            print 
            
            for o in r4:
                print r4[o] 
                 
        
            self.assertEqual(len(r3),2)

        
if __name__ == '__main__':
    BT.localTest()    
