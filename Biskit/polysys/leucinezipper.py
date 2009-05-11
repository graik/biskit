import Biskit.molUtils as MU
from numpy import correlate,array, corrcoef,where

class LeuZip:
    props = [   ["LEU","ILE","ALA"], ##a , nonpolar
                ["ILE","PHE","VAL","LEU","TRP","MET","ALA","GLY","CYS","TYR","PRO"], ##b , hidrophilic
                ["ILE","PHE","VAL","LEU","TRP","MET","ALA","GLY","CYS","TYR","PRO"], ##c , hidrophilic
                ["LEU","ILE","ALA"], ##d , nonpolar
                ["GLU","GLN"], ##e , polar electrostatic
                ["ILE","PHE","VAL","LEU","TRP","MET","ALA","GLY","CYS","TYR","PRO"], ##f , hidrophilic
                ["GLU","GLN","ARG","LYS"]] ##g , polar electrostatic
    
    scores = [  [1.5,.9,.8], ##a , nonpolar
                [.8,.75,.7,.65,.6,.55,.5,.4,.3,.2,.1], ##b , hidrophilic
                [.8,.75,.7,.65,.6,.55,.5,.4,.3,.2,.1], ##c , hidrophilic
                [1.,.9,.8], ##d , nonpolar
                [.9,.9], ##e , polar electrostatic
                [.8,.75,.7,.65,.6,.55,.5,.4,.3,.2,.1], ##f , hidrophilic
                [.8,.8,.9,.9]] ##g , polar electrostatic
    
    window_length = 7
    
    def __init__(self):
        pass
    

    
    def findHeptads(self, chain, k=20):
        tries  = 0
        myscores = self.copyScores()
        heptads = [""]*k
        corr = []*k
        while tries < k:
            score = self.scoreHeptads(chain,myscores)
            heptads[tries] = self.heptadize(score,chain)
            #~ print heptads[tries]
            myscores = self.lowContributions(heptads[tries],myscores)
            lol = 0
            for m in score:
                print lol,m
                lol = lol+1
            #~ corr[tries] = correlate(chain,heptads[tries],mode= "full")
            tries = tries+1
    
    def heptadize(self,score,chain):
        s = array(score)
        top = max(score)
        index = where(s == top)[0][0]
        #~ print top
        if index +7 >= len(chain):
            return ""
        return  chain[index:index+7]
    
    def scoreHeptads(self, chain, scores):
        window = 0
        #~ print chain, len(chain)
        score = [0]*len(chain)
        while window + 7 <= len(chain):
            score[window] = 0
            for i in range(0,self.window_length):  
                for j in range(0,len(self.props[i])):
                    if MU.single2longAA(chain[window+i])[0] == self.props[i][j]:
                        score[window] = score[window] + scores[i][j]
            window = window+1
            
        #~ for m in score:
            #~ print m
        #~ print score
        #~ print len(correlate(score,score,mode= "full"))
        #~ for m in correlate(score,score,mode= "full"):
            #~ print m
        return score
        
        
    def copyScores(self):
        copy = []
        for i in self.scores:
            copy.append([])
            for j in i:
                copy[-1].append(j)
        #~ print copy
        return copy
        
        
    def lowContributions(self,heptad, scores):
        window = 0
        for c in heptad:
            for k in range(0,len(self.props[window])):
                #~ print MU.single2longAA(c) , self.props[window][k], window
                if  MU.single2longAA(c)[0] == self.props[window][k]:
                    scores[window][k] =  scores[window][k] - .1
                    #~ print "*"
            window = window +1
        #~ print scores
        return scores

##############
## Test
##############
import Biskit.test as BT
import Biskit.tools as T
from Biskit.PDBModel import PDBModel
import os

class Test(BT.BiskitTest):
    """ Test cases for Polyfret"""

    def prepare(self):
        pass

    def cleanUp( self ):
        pass
        
    def test_search(self):
        """LeuZip pattern search test cases"""

        l = LeuZip()
        l.findHeptads("LEIRAAFLRRRNTALRTRVAELRQRVQRLRNIVSQYETRYGPL")
    


if __name__ == '__main__':
    BT.localTest()    