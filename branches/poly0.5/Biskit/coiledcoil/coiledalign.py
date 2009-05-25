from coiledcoil import CoiledCoil
import Biskit.molUtils as MU
from coiledutils import scores2String, scoreAlignFun
from Bio.pairwise2 import align
from numpy import array
import Biskit.molUtils as MU


class CoiledAlign:
    """
    Class for coiled coil parameters study. Mainly it will be used for getting the better alignment
    for an unknown cc structure to do homology modelling.
    """
    
    
    
    def __init__(self,cc = None, likeness_table = "BLOSSUM62"):
        
        self.like_scores, self.charge_scores = self.parseLikelyhood(likeness_table)
        self.mycc = cc or CoiledCoil()

    
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
    
    
    
    def alignChains(self, a = "", b = "", reg_a = "", reg_b = ""):
        """
        
        """
        assert(len(a) == len(reg_a) and len(b) == len(reg_b) ),"No correspondence between chains and their registers."
        
       
        ## Score alignment
        scores_a = []
        scores_b = []
        
        pos = "abcdefg".find(reg_a[0])
        for c in a:
            scores_a.append(self.mycc.scores[MU.single2longAA(c)[0]][pos])
            pos = (pos+1) % 7
        
        pos = "abcdefg".find(reg_b[0])
        for c in b:
            scores_b.append(self.mycc.scores[MU.single2longAA(c)[0]][pos])
            pos = (pos+1) % 7
       
        a_s,b_s = scores2String(scores_a,scores_b)

        scc = align.localcs(a_s,b_s,scoreAlignFun,-1000,-1000)[0][2:4]
        #~ print align.localcs(a_s,b_s,scoreAlignFun,-1000,-1000)
        
        ## Chain similarity
       
        chs = align.localds(a, b, self.like_scores,-1000,-1000)[0][2:4]
        
        ## Charges similarity
        crs = align.localds(a, b, self.charge_scores,-1000,-1000)[0][2:4]
        
        return {"heptad":(self.getFitnessScore("heptad",a,b,scc[1],reg_a,reg_b),scc[1]) \
            , "res_like":(self.getFitnessScore("res_like",a,b,chs[1],reg_a,reg_b),chs[1]) \
            , "charges": (self.getFitnessScore("charges",a,b,crs[1],reg_a,reg_b),crs[1])}
    
    def getFitnessScore(self,type = "heptad",a = "",b = "",where = 0,reg_a = "",reg_b=""):
        """
        "heptad"
        "res_like"
        "charges"
        """
        
        assert( type in ["heptad","res_like","charges"] ), "Please choose a correct type (heptad,res_like,charges)." 
        
        chain_a = a
        chain_b = b
        freg_a = reg_a
        freg_b = reg_b
        
        if type ==  "heptad":
            assert(len(a) == len(reg_a) and len(b) == len(reg_b) ),"No correspondence between chains and their registers."
        
        
        a = a[where:where+len(b)]
        reg_a = reg_a[where:where+len(b)]
        
        if len(a)< len(b):
            b = b[0:len(a)]
            reg_b = reg_b[0:len(a)]
    
        if type ==  "heptad":
            scores_a = []
            scores_b = []
            
            pos = "abcdefg".find(reg_a[0])
            for c in a:
                scores_a.append(self.mycc.scores[MU.single2longAA(c)[0]][pos])
                pos = (pos+1) % 7
            
            pos = "abcdefg".find(reg_b[0])
            for c in b:
                scores_b.append(self.mycc.scores[MU.single2longAA(c)[0]][pos])
                pos = (pos+1) % 7
       
            #~ ## get scores
            #~ scores_af = self.flatten(chain_a,freg_a)
            #~ scores_bf = self.flatten(chain_b,freg_b)
            #~ scores_a = self.flatten(a,reg_a)
            #~ scores_b = self.flatten(b,reg_b)
            
            #~ chain_as ,chain_bs = scores2String(scores_a+scores_af,scores_b+scores_bf) 
            
            chain_as ,chain_bs = scores2String(scores_a,scores_b)
            
            #~ chain_as2 = chain_as[:len(scores_a)]
            #~ chain_bs2 = chain_bs[:len(scores_b)]
            
            total_score = 0
            for i in range(len(chain_as)):
                total_score += scoreAlignFun(chain_as[i],chain_bs[i])
            
            return total_score
            
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
        
        start = register.find("a")
        end  = start + len(chain[start:]) - len(chain[start:])%window_length
      
        
        scores = []
        pos = 0
        ac = 0 

        for i in chain[start:end]:
            ac+= self.mycc.scores[MU.single2longAA(i)[0]][pos] or 0
            pos = (pos+1)%window_length
            if pos == 0:
                scores.append(ac)
                ac = 0
        
        return scores
    
    def alignRegisters(self,a = "",b = "",reg_a ="",reg_b =""):
        """
        """
        ch_a=a;ch_b=b;bra=reg_a;brb=reg_b
        
        start_a = reg_a.find("a")
        start_b = reg_b.find("a")
        ## So we make it end with reg letter g
        a = a[start_a:]
        
        if len(a)%7 != 0:
            a = a[:(len(a)/7)*7]
        
        b = b[start_b:]
        if len(b)%7 != 0:
            b = b[:(len(b)/7)*7]
       
        
        total_regs_a = len(a) / 7
        total_regs_b = len(b)  / 7
        
        ## Scores
        scores_a = self.flatten(a,"abcdefg"*total_regs_a)
        scores_b = self.flatten(b,"abcdefg"*total_regs_b)
        
        a_s,b_s = scores2String(scores_a,scores_b)
        
        scc = align.localcs(a_s,b_s,scoreAlignFun,-1000,-1000)[0][2:4]

        scc =( scc[0],start_a+scc[1]*7)
        
        
        
        ## Charges and residues
        subvals_res = []
        subvals_charge = []
        for i in range(total_regs_a - total_regs_b +1):
            pos = i*7
            acc_charges = 0
            acc_res = 0
            for j in range(len(b)):
                acc_res += self.like_scores[(a[pos+j],b[j])]
                acc_charges += self.charge_scores[(a[pos+j],b[j])]
            subvals_charge.append((acc_charges,start_a+pos))
            subvals_res.append((acc_res,start_a+pos))
        
        crs = max(subvals_charge)
        chs = max(subvals_res)
        
        #~ return {"heptad": scc, "charges": max(subvals_charge),"res_like": max(subvals_res)}
        a = ch_a;b=ch_b;reg_a=bra;reg_b=brb
        return {"heptad":(self.getFitnessScore("heptad",a,b,scc[1],reg_a,reg_b),scc[1]) \
            , "res_like":(self.getFitnessScore("res_like",a,b,chs[1],reg_a,reg_b),chs[1]) \
            , "charges": (self.getFitnessScore("charges",a,b,crs[1],reg_a,reg_b),crs[1])}
    
        
##############
## Test
##############
import Biskit.test as BT
import Biskit.tools as T


class Test(BT.BiskitTest):
    """ Test cases for Coiled Coil Alignment"""

    def prepare(self):
        self.cc = CoiledCoil(T.dataRoot() + '/coiledcoil/SOCKET_par_norm')
        self.ca = CoiledAlign(self.cc)

    def cleanUp( self ):
        pass
        
    def test_align(self):
        """Best alignment search test cases"""
        res = self.ca.alignChains("RRRLLLLLLLRRRRRRLLLLLLLRRR","LLLRLRL","efgabcdefgabcdefgabcdefgab","abcdefg")
        self.assertEqual(res,{'res_like': (18.0, 7), 'charges': (3.6000000000000001, 20), 'heptad': (69898, 7)})
    
    
    def test_align_regs(self):
        """Best alignment for registers test cases"""
        
        res = self.ca.alignRegisters("RRRLLLLLLLRRRRRRLLLLLLLRRR", "LLLRRRR","efgabcdefgabcdefgabcdefgab","abcdefg")
        self.assertEqual( res,{'res_like': (11.0, 17), 'charges': (2, 10), 'heptad': (69890, 3)})
        res = self.ca.alignRegisters("RRRLLLLLLLRRRRRRLLLLLLLRRR", "RRRLLLL","efgabcdefgabcdefgabcdefgab","abcdefg")
        self.assertEqual(res,{'res_like': (13.0, 10), 'charges': (3.2000000000000002, 10), 'heptad': (69975, 3)})
        res = self.ca.alignRegisters("RRRLLLLLLLRRRRRRLLLLLLLRRR", "RRRRRRR","efgabcdefgabcdefgabcdefgab","abcdefg")
        self.assertEqual( res,{'res_like': (28.0, 10), 'charges': (11, 10), 'heptad': (69865, 3)})
        
    def test_fitness(self):
        """getFitnessScore function test cases"""
        self.assertEqual( self.ca.getFitnessScore("heptad","RRRLLLLLLLRRRRRRLLLLLLLRRR","LLLRLRL",7,"efgabcdefgabcdefgabcdefgab","abcdefg"),69898)
        self.assertEqual( self.ca.getFitnessScore("res_like","RRRLLLLLLLRRRRRRLLLLLLLRRR","LLLRLRL",7,"efgabcdefgabcdefgabcdefgab","abcdefg"),18.0)
        self.assertEqual( self.ca.getFitnessScore("charges","RRRLLLLLLLRRRRRRLLLLLLLRRR","LLLRLRL",20,"efgabcdefgabcdefgabcdefgab","abcdefg"),3.6)
        
       
if __name__ == '__main__':
    BT.localTest()    
    