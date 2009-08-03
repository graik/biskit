##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2009 Victor Gil & Raik Gruenberg
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.
##
## You find a copy of the GNU General Public License in the file
## license.txt along with this program; if not, write to the Free
## Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
##
##
## last $Author$
## last $Date: 2009-03-09 09:10:46 +0100 (Mon, 09 Mar 2009) $
## $Revision$

import Biskit.molUtils as MU
from Biskit import BiskitError
from numpy import correlate,array, where, sum
from coiledutils import getRegister,flatten



class CoiledCoil:
    """
    Temptative class for Coiled Coil registry analysis.
    An heptad registry is given by a sequence of letters from a to g,
    where a corresponds to the first residue and g to the last.
    Ex.
    
    LQAIEKQ
    abcdefg
    """
    
    
    
    """
    Size of the window. While it's generally 7, another repeating structures
    can happen so this is left as a constant just for completeness sake.
    (Not used)
    """
    window_length = 7 
    
    
    
    def __init__(self,db = "",find_style = "mean"):
        """
        Instantiation function. It loads a score table from disk. Lines of the score 
        table must have this format:
        
        (aminoacid 3 letter code) (probability of this aa to be in position 'a') (prob of being in position b)...
        
        Ex.
        
        LEU 0.321 0.035 0.045 0.311 0.067 0.034 0.045
        
        @param db: Path of the scores file.
        @type db: string
        @param find_style: Method used for selecting the best heptad. It can be:
            - "mean" to use the mean of all scores.
            - "max" to use the maximun score.
        @type find_style: string
        
        """
        self.db = db or T.dataRoot() + '/coiledcoil/DADParry_scaled'

        self.scores = self.parseScores( self.db )
        
        self.find_style = "mean"
        
        
        
    def parseScores( self, db="" ):
        """
        Load position scores from external file (probability for an aa X to be
        in register position Y, so it's a 20x7 table).
        
        @param db: path to file with amino acid scores per position
        @type  db: str
        """
        r = {}
        db = db or self.db
        ## Load data from file
        try:
            lineas = open(db,"r").readlines()
        except IOError, msg:
            raise BiskitError('cannot open score file %s.\n Error: %s' \
                              % (db, msg ) )
        
        lineas = [ l.strip() for l in lineas ]
        
        for l in lineas:
            aux = l.split()
            if len(aux[0]) == 3 :
                r[aux[0]] = [ float( n ) for n in aux[1:] ]
            else:
                r[MU.single2longAA(aux[0])[0]] = [ float( n ) for n in aux[1:] ]
            
        return r

    
    def findHeptads(self, chain = "", k=10):
        """
        Function for discovering the heptad register of a coiled coil protein.
        
        @param chain: Chain to be analyzed with aminoacids in single-letter code.
        @type chain: string
        @param k: Number of heptads to try to fit into the chain (if 1 it will only
                try to fit the best scored one, if 2 it will also try with the second 
                best scored... and so on. Best scored heptad is not allways the target
                heptad, so k > 1 is recommended.
        @type k: integer
        
        @return: A dictionary with entries:
                - "heptad" for the best heptad (a tuple with score (int) and heptad chain (string)).
                - "reg" for the registry sequence for the whole chain.
                - "corr" for the accumulated correlation for the best scored heptad.
        @rtype: Dictionary( "best":tuple(float,string), "reg":string, "corr":list(float))
        
        """
        #~ print self.find_style, self.db
        
        tries  = 0
        
        score = self.scoreHeptads(chain)
               
        heptads = self.heptadize(score,chain,k)
        
        if self.find_style == "max":
            indexes = {}
            new_scores = [[0.,0.]]*len(heptads)
            for i in range(len(heptads)):
                index = heptads[i][1]
                indexes[heptads[i][0]] = []
                new_scores[i] = [0.,heptads[i][0]]
                while index <len(chain):
                    new_scores[i][0] += score[index]
                    indexes[heptads[i][0]].append(index) 
                    index += self.window_length
            
            for i in range(len(heptads)):
                index = heptads[i][1]- self.window_length
                while index >=0:
                    new_scores[i][0] += score[index]
                    indexes[heptads[i][0]].append(index) 
                    index -= self.window_length
                indexes[heptads[i][0]].sort()
           
            self.last_heptads = []
            for i in range(len(new_scores)):
                self.last_heptads.append( (new_scores[i][0],heptads[i]))
            self.last_score = score
            
            
            ## And the winner is...
            best = sorted(new_scores,reverse=True)[0][1]
            self.indexes = indexes[best]
        else:
            #~ print "*"
            best = sorted(score,reverse=True)[0][1]
        
        c = self.correlate(chain,[best])[1]
        
        
        
        retdic = {}
        retdic["best"] = best
        retdic["reg"] = getRegister(best,chain)
        retdic["corr"] = c
        
        return retdic
    
    
        
    def correlate(self, chain, heptads):
        """
        Function for calculating the accumulated correlation of a chain for a 
        given set of heptads.
        
        @param chain: Chain to be correlated with the heptads.
        @type chain: string
        @param heptads: List containing the heptads we want to get the correlation
            with @chain.
        @type heptads: List (string)
        
        
        @return: It returns a tuple of two structures:
            
            - Corr: Which is a dictionary indexed by the string representing
                the residues of the heptad (so these inside @heptads) and containing
                for each position the peak of correlation
            
            while...
            
            - Acccum: Is another dictionary indexed in the same way containing the 
                accumulated correlation values for each position. Usefull to see 
                repeating structures.
        
        @rtype: Tuple(Dictionary,Dictionary)
        """
        
        corr = {}
        accum = {}
        for k in range(0,len(heptads)):
            aux= []
            corr[heptads[k]] = []
            for j in range(7):
                (a,b,c) = self.sequenceCorrelation(chain,j,heptads[k])
                corr[heptads[k]].append(b)
                if aux == []:
                    aux = c
                else:
                    aux+=c
            accum[heptads[k]] = aux
        
        return corr , accum
    
    
    def sequenceCorrelation	(self, a="",channel=0,b = ""):
        """
        Version of polysys::Protein::sequenceCorrelation. It calculates the 
        score correlation between 2 peptide chains.
        
        @param a: Chain of aminoacids in single-letter code.
        @type a: string
        @param channel: Registry frame of chain a to be used.
        @type channel: integer
        @param b: Chain of aminoacids in single-letter code.
        @type b: string
        
        @return: It returns the peak of correlation and the index where this peak
                happens. Last value returned is the correlation list.
        @rtype: tuple (float, [integer], list(float) )
        """
        
        al = []
        bl = []
        
        pos= channel
        for i in a:
            al+= [self.scores[MU.single2longAA(i)[0]][pos] or 0]
            pos = (pos+1)%self.window_length
            
        pos = 0
        for i in b[0]:
            bl+= [self.scores[MU.single2longAA(i)[0]][pos] or 0]
            pos = pos+1

        a1 = array(al)
        b1 = array(bl)
        ab = correlate(a1,b1)
        
        
        return (max(ab) , where(ab==max(ab)),ab)	
    
    
    def heptadize(self,score,chain,k):
        """
        Function for retrieving the best k heptads from score calculations.
        
        @return: List of the k best heptads.
        @rtype: list
        """
        
        l = []
        for i in range(0,len(score)):
            l.append((score[i],i))
        
        s = sorted(l,reverse=True)
        heptads = []
        for e in range(0,k):
            index = s[e][1]
            if not (index +7 >= len(chain)):
                heptads.append((chain[index:index+7],index))
        return heptads
            
        
            
    def scoreHeptads(self, chain):
        """
        Function to calculate the score of each window in the chain.
        
        @return: List of scores for each position.
        @rtype: list
        
        """
        if self.find_style == "max":
            window = 0
            score = [0.]*len(chain)
            while window + 7 <= len(chain):
                score[window] = 0
                for i in range(0,self.window_length):  
                    r = MU.single2longAA(chain[window+i])[0]
                    if  r in self.scores:
                        score[window]+=self.scores[r][i]
                window = window+1
            #~ print score
            return score
        
        elif self.find_style == "mean":
            print "*",
            candidates = []
            score = []
            
            for i in range(0,len(chain)):
                if i+7 <=len(chain):
                    candidates.append(chain[i:i+7])
            #~ print chain
            #~ print candidates
            for c in candidates:
                reg = getRegister(c,chain)
                scores = flatten(chain,reg,self)
                total = 0.
                for s in scores:
                    total += s
                mean = total /len(scores)
                score.append((mean,c))
            
            score += [0.,""]*(len(chain)-len(score))
            #~ print score, len(score),len(chain),len(candidates)
            return score
        else:
            print "Error in CoiledCoil:%s is not an allowed \"find style\" method"%self.find_style 
            return []

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
        """Pattern search test cases"""
        l = CoiledCoil()
        l2 = CoiledCoil(T.dataRoot() + '/coiledcoil/SOCKET_par_norm')
        self.assertEqual(len(l.scores),20)
        if  l.find_style == "max":
            
            
            ##TARGET
            target = 'LEIRAAFLRRRNTALRTRVAELRQRVQRLRNIVSQYETRYGPL'
            res = l.findHeptads(target)
            
            if self.local :
                print target
                print res["reg"]
                for i in l.last_heptads:
                    print i[0],i[1]
            
            
            self.assertEqual(res["best"],"LRTRVAE")
            self.assertEqual(len(res["reg"]),43)
            
            res = l2.findHeptads(target)
            
            if self.local:
                print target
                print res["reg"]
                for i in target :
                    print i,"   ",
                print
                
                for i in l2.last_score:
                    print "%.2f " %(i),
                print
            
            self.assertEqual(res["best"],"VSQYETR")
            self.assertEqual(len(res["reg"]),43)

    def test_correlation(self):
        """Correlation test case"""
        print "WARNING!!!!!!!!! & REMEMBER!!!!! You HAVE TO implement this test case, as correlation isn't fully tested yet!!!!"
    
    def test_mean_scoring(self):
        """Mean scoring vs max scoring test case"""
        l = CoiledCoil()
        print
        chain = ""
        for i in range(len(MU.allAA())):
            chain += MU.allAA()[i]
        print chain
        l.scoreHeptads(chain+chain)
        print l.findHeptads(chain+chain)
        l.find_style = "max"
        l.scoreHeptads(chain+chain)

if __name__ == '__main__':
    BT.localTest()    
    
    
    
    
