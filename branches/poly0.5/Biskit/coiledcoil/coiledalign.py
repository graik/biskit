


class CoiledAlign:
    """
    Class for coiled coil parameters study. Mainly it will be used for getting the better alignment
    for an unknown cc structure to do homology modelling.
    """
    
    
    
    def __init__(self,likeness_table = "BLOSSUM62"):
        
        self.like_scores, self.charge_scores = self.parseLikelyhood(likeness_table)
    

    
    def parseLikelyhood(self,table = ""):
        """
        Load residue likelyhood scores for BLOSSUM 62 table.
        """
        try:
            lineas = open(T.dataRoot() + '/coiledcoil/'+table,"r").readlines()
        except IOError, msg:
            raise BiskitError('cannot open score file %s.\n Error: %s' \
                              % (db, msg ) )
        
        like_scores = {}
        charge_scores = {}
        
        lineas = [ l.strip() for l in lineas ]
        
        residues = "A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *".split()
        
        charged_residues = "E R H D E".split()
        
        for l in lineas:
            aux = l.split()
            scores = [ float( n ) for n in aux[1:] ]
            for i in range(len(scores)):
                like_scores[(aux[0],residues[i])] = scores[i]
                if aux[0] in charged_residues or residues[i] in charged_residues :
                    charge_scores[(aux[0],residues[i])] = -1
                else:
                    charge_scores[(aux[0],residues[i])] = 0.2
        
        charge_scores[("K","K")] = 2
        charge_scores[("R","R")] = 2
        charge_scores[("H","H")] = 2
        charge_scores[("D","D")] = 2
        charge_scores[("E","E")] = 2
        
        charge_scores[("K","R")] = charge_scores[("R","K")] = 1.5
        charge_scores[("K","H")] = charge_scores[("H","K")] = 1.5
        charge_scores[("R","H")] = charge_scores[("H","R")] = 1.5
        charge_scores[("D","E")] = charge_scores[("E","D")] = 1.5
        
        charge_scores[("D","K")] = charge_scores[("K","D")] = -2
        charge_scores[("E","K")] = charge_scores[("K","E")] = -2
        charge_scores[("D","R")] = charge_scores[("R","D")] = -2
        charge_scores[("E","R")] = charge_scores[("R","E")] = -2
        charge_scores[("D","H")] = charge_scores[("H","D")] = -2
        charge_scores[("E","H")] = charge_scores[("H","E")] = -2
        
        return like_scores, charge_scores
    
    
    
    def tryAlignment(self, a = "", b = "", reg_a = "", reg_b = ""):
        """
        
        """
        if len(a) < len(b):
            aux =  b
            b  = a
            a = aux
            aux =  reg_b
            reg_b  = reg_a
            reg_a = aux
        
        ## Score correlation
        scores_a = self.flatten(a,reg_a)
        scores_b = self.flatten(b,reg_b)
       
        a,b = self.toStringChain(scores_a,scores_b)

        scc = align.localcs(a,b,self.scoreAlignFun,-1000,-1000)[0][2:4]
        
        ## Chain similarity
        chs = align.localds(a, b, self.like_scores,-100,-100)[0][2:4]

        ## Charges similarity
        crs = align.localds(a, b, self.charge_scores,-100,-100)[0][2:4]
        
        return {"heptad":scc , "res_like":chs , "charges":crs}
    
    def getFitnessScore(self,type = "heptad",a = "",b = "",where = 0,reg_a = "",reg_b=""):
        
        """
        "heptad"
        "res_like"
        "charges"
        """
        
        assert( type in ["heptad","res_like","charges"] ), "Please choose a correct type (heptad,res_like,charges)." 
        print len(a) ,len(reg_a),  len(b), len(reg_b)
        
        if type ==  "heptad":
            assert(len(a) == len(reg_a) and len(b) == len(reg_b) ),"No correspondence between chains and their registers."
        
        a = a[where:where+len(b)]
        reg_a = reg_b[where:where+len(b)]
        
        if len(a)< len(b):
            b = b[0:len(a)]
            reg_b = reg_b[0:len(a)]
    
        if type ==  "heptad":
        
            ## get scores
            a_s = []
            b_s = []
            
            pos = reg_a.find("abcdefg")
            for c in a:
                a_s.append(self.scores[MU.single2longAA(c)[0]][pos])
                pos = (pos+1)%self.window_length
            
            
            pos = reg_b.find("abcdefg")
            for c in b:
                b_s.append(self.scores[MU.single2longAA(c)[0]][pos])
                pos = (pos+1)%self.window_length
            
            a_sa = array(a_s)
            b_sb = array(b_s)
            
            return (a_sa*b_sb).sum() 
            
        elif type == "res_like":
            acc = 0
            for i in range(len(a)):
                acc += self.like_scores[(a[i],b[i])]
            return acc
            
        else:
            acc = 0
            for i in range(len(a)):
                acc += self.charge_scores[(a[i],b[i])]
            return acc
            
            
    def flatten (self, chain = "", register = "", window_length = 7):
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
        #~ print chain , register
        start = register.find("a")
        end  = start + len(chain[start:]) - len(chain[start:])%window_length

        scores = []
        pos = 0
        ac = 1 ## So we don't want values smaller than 1 

        for i in chain[start:end]:
            ac+= self.scores[MU.single2longAA(i)[0]][pos] or 0
            pos = (pos+1)%self.window_length
            if pos == 0:
                scores.append(ac)
                ac = 1

        return scores
        
##############
## Test
##############
import Biskit.test as BT
import Biskit.tools as T
from coiledcoil import CoiledCoil

class Test(BT.BiskitTest):
    """ Test cases for Coiled Coil Alignment"""

    def prepare(self):
        self.ca = CoiledAlign()

    def cleanUp( self ):
        pass
        
    def test_align(self):
        """Best alignment search test cases"""
        res = self.ca.tryAlignment("RRRLLLLLLLRRR","LRLRLRL","efgabcdefgabc","abcdefg")
        print res
        
         
    def test_fitness(self):
        """getFitnessScore function test cases"""
        #~ self.assertEqual( l2.getFitnessScore("heptad","XXXLLLLLLXXX","AAAAAA",3,"cdefgab"),7.0316)
        self.assertEqual( self.ca.getFitnessScore("res_like","XXXLLLLLLLXXX","AAAAAAA",3),-7.0)
        self.assertEqual( self.ca.getFitnessScore("charges","XXXLLLLLLLXXX","AAAAAAA",3),1.4)
        
       
if __name__ == '__main__':
    BT.localTest()    
    