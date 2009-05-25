
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
        print
        
        for h in heptads.keys():
            score = {}
            #~ print h
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
            
            #~ print score    
            priorities = {"CC":4,"Pa":2,"Pair":4,"Paper":5,"So":3}

            best = (score.keys()[0],0)
            for k2 in score.keys():
                if best[1] < score[k2]:
                    best = (k2,score[k2])
                if best[1] == score[k2]: ## In this case we have to manage priorities
                    if priorities[k2] > priorities[best[0]]:
                        best  = (k2,score[k2])
            
            heptads[h]["best"] = best
            #~ print best
        
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
            
        scores = {}
        for k in hits.keys():
            scores[k] = float(hits[k])/tries[k]
        
        
        print scores
        
        
        ## Get Alignment values
        
        
        target_hep = heptads["TARGET"][heptads["TARGET"]["best"][0]]
        if heptads["TARGET"]["best"][0]== "Paper":
            target_hep = heptads["TARGET"][heptads["TARGET"]["best"][0]][0]
        target = heptads["TARGET"]["seq"]
        target_reg=getRegister(target_hep,target)
        
        
        cross_scores_chain = {}
        cross_scores_reg = {}
        
        for h in ["2Z5H"]:##heptads.keys():
            print h
            other_hep = heptads[h][heptads[h]["best"][0]]
            if heptads[h]["best"][0]== "Paper":
                other_hep = heptads[h][heptads[h]["best"][0]][0]
            other = heptads[h]["seq"]
            other_reg=getRegister(other_hep,other)
            #~ print "l",len(other_reg),len(other)
            #~ print other, target, other_reg, target_reg
            
            res = ca.alignChains(other, target, other_reg, target_reg)
            res2 = ca.alignRegisters(other, target, other_reg, target_reg)
            print res
            print res2
            
            cross_scores_chain[h]={}
            for i in ["heptad","res_like","charges"]:
                for j in ["heptad","res_like","charges"]:
                    cross_scores_chain[h][(i,j)] = ca.getFitnessScore(j,other,target,res[j][1], other_reg,target_reg)
                
            print cross_scores_chain[h]
            
            cross_scores_reg[h]={}
            for i in ["heptad","res_like","charges"]:
                for j in ["heptad","res_like","charges"]:
                    cross_scores_reg[h][(i,j)]=  ca.getFitnessScore(j,other,target,res2[j][1], other_reg,target_reg)
                
            print cross_scores_reg[h]
        
        
        ## primero eliminar los offsets de 0, luego normalizar a 1
        
        #~ for k in ["2Z5H"]:##cross_scores_chain:
            #~ for i in ["heptad","res_like","charges"]:
                #~ mymax = cross_scores_chain[k][(i,"heptad")]
                #~ for j in ["heptad","res_like","charges"]:
                    #~ if mymax < cross_scores_chain[k][(i,j)]:
                        #~ mymax = cross_scores_chain[k][(i,j)]
                #~ for j in ["heptad","res_like","charges"]:
                    #~ cross_scores_chain[k][(i,j)] /= mymax

                
       
        for k in ["2Z5H"]:##cross_scores_chain:
            print "\n        h\tr\tc"
            for i in ["heptad","res_like","charges"]:
                print "%-8s"%(i),
                for j in ["heptad","res_like","charges"]:
                    print "%.2f / %.2f"%(cross_scores_chain[k][(i,j)],cross_scores_reg[h][(i,j)] ),
                print
            
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
        
if __name__ == '__main__':
    BT.localTest()    
       