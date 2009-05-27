
from coiledcoil import CoiledCoil
from coiledalign import CoiledAlign
from coiledutils import sameRegister, getRegister

class CCStudy:
    
    def __init__(self,data ="",cc = None):
        self.parseData(data)
        self.mycc = cc or CoiledCoil() 
        
    def parseData(self,dat = ""):
        
        try:
            lineas = open(dat,"r").readlines()
        except IOError, msg:
            raise BiskitError('cannot open score file %s.\n Error: %s' \
                              % (db, msg ) )
        
        lineas = [ l.strip() for l in lineas ]
        
        self.data = {}
        
        for l in lineas:
            
            aux = l.split()
            if aux!=[]:
                self.data[aux[0]]= {}
                for w in aux[1:]:
                    aux2 = w.split(":")
                    if aux2[0]!="Paper":
                        self.data[aux[0]][aux2[0]]= aux2[1]
                    else:
                        self.data[aux[0]][aux2[0]]= (aux2[1],aux2[2])
            
    
    def doStudy(self):
        ca = CoiledAlign(self.mycc)
        
        
        ## Statistics
        sources = ["CC","Pa","Pair","Paper","So"]
        hits = {}
        tries = {}
        fails = {}
        for k in sources:
            hits[k] = 0
            tries[k] = 0
        sources.append("seq")
        
        heptads = {}
        for d in self.data.keys():
            ## Extract heptads
            heptads[d]={}
            for k in self.data[d].keys():
                if k in sources:
                    heptads[d][k] = self.data[d][k]
        
        
        for h in heptads.keys():
            score = {}
            there = heptads[h].keys()
            there.remove("seq")
            for k in there:
                com_keys = heptads[h].keys()
                com_keys.remove(k)
                com_keys.remove("seq")
                score[k] = 0
                for k2 in com_keys:
                    if k != "Paper":
                        one = heptads[h][k]
                    else:
                        one = heptads[h][k][0]
                    if k2 != "Paper":
                        the_other = heptads[h][k2]
                    else:
                        the_other = heptads[h][k2][0]
                    
                    if sameRegister(one,the_other,heptads[h]["seq"]):
                        score[k] += 1
            
            
            priorities = {"CC":4,"Pa":2,"Pair":4,"Paper":5,"So":3}

            best = (score.keys()[0],0)
            for k2 in score.keys():
                if best[1] < score[k2]:
                    best = (k2,score[k2])
                if best[1] == score[k2]: ## In this case we have to manage priorities
                    if priorities[k2] > priorities[best[0]]:
                        best  = (k2,score[k2])
            
            heptads[h]["best"] = best
            
        
        for h in heptads.keys():
            there = heptads[h].keys()
            there.remove("seq")
            there.remove("best")
            for k in there:
                tries[k] +=1
                
                one = heptads[h][k]
                if k == "Paper":
                    one = heptads[h][k][0]
                    
                the_other = heptads[h][heptads[h]["best"][0]]
                if heptads[h]["best"][0] == "Paper":
                    the_other = heptads[h][heptads[h]["best"][0]][0]
                
                if sameRegister(one,the_other,heptads[h]["seq"]):
                    hits[k] += 1
            
        self.scores = {}
        for k in hits.keys():
            self.scores[k] = (float(hits[k])/tries[k],hits[k],tries[k])
        
        
        print self.scores
        
        
        ## Get Alignment values
        
        
        target_hep = heptads["TARGET"][heptads["TARGET"]["best"][0]]
        if heptads["TARGET"]["best"][0]== "Paper":
            target_hep = heptads["TARGET"][heptads["TARGET"]["best"][0]][0]
        target = heptads["TARGET"]["seq"]
        target_reg=getRegister(target_hep,target)
        
        
        self.alignments = {}
        
        for h in ["DUMMY"]:##heptads.keys():
            ## Vigilar aqui cual escoger segun la tabla!!
            self.alignments[h] = ca.copy()
            print h
            other_hep = heptads[h][heptads[h]["best"][0]]
            if heptads[h]["best"][0]== "Paper":
                other_hep = heptads[h][heptads[h]["best"][0]][0]
            other = heptads[h]["seq"]
            other_reg=getRegister(other_hep,other)
                
            self.alignments[h].alignChains(other, target, other_reg, target_reg)
            self.alignments[h].alignRegisters(other, target, other_reg, target_reg)
            
        return self.scores, self.alignments
            
    
    def chooseBest(self):
        ## First for each choose the best
        maxims = []
        for k in self.alignments.keys():
            print "key: ",k
            self.alignments[k].normalizeScores()
            keys = self.alignments[k].chain_alignments.keys()
            maxacc = (0,0)
            print keys[0]
            for t in self.alignments[k].chain_alignments[keys[0]]:
                print t[0],t[1]
                acc = ((t[0] + self._findInAlignment(self.alignments[k].chain_alignments[keys[1]],t[1])+ self._findInAlignment(self.alignments[k].chain_alignments[keys[2]],t[1]),t[1]))
                maxacc = max(maxacc,acc)
            maxims.append((maxacc,k))
        ## Then choose the best of all !!!
        best_chain = max(maxims)
        
        maxims = []
        for k in self.alignments.keys():
            keys = self.alignments[k].reg_alignments.keys()
            maxacc = (0,0)
            for t in self.alignments[k].reg_alignments[keys[0]]:
                print t[0],t[1]
                acc = ((t[0] + self._findInAlignment(self.alignments[k].reg_alignments[keys[1]],t[1])+ self._findInAlignment(self.alignments[k].reg_alignments[keys[2]],t[1]),t[1]))
                maxacc = max(maxacc,acc)
            maxims.append((maxacc,k))
        ## Then choose the best of all !!!
        best_reg = max(maxims)
        
        return best_chain, best_reg
    
    def _findInAlignment(self,where=[],what = 0):
        for i in range(len(where)):
            if where[i][1] == what:
                return where[i][0]
        return 0
        
    def createDataFile(self, struct_file = "",file="a.dat"):
        
        pass
        

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
        
    def test_Study(self):
        """doStudy function test case"""
        cs = CCStudy(T.dataRoot() + '/coiledcoil/coils.dat')
        cs.doStudy()
        
        print cs.chooseBest()
if __name__ == '__main__':
    BT.localTest()    
       