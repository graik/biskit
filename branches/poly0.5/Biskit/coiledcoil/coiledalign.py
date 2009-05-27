from coiledcoil import CoiledCoil
from coiledutils import scores2String, scoreAlignFun,alignd,alignf,flatten
import Biskit.molUtils as MU



class CoiledAlign:
    """
    Class for coiled coil parameters study. Mainly it will be used for getting the better alignment
    for an unknown cc structure to do homology modelling.
    """
    
    
    
    def __init__(self,cc = None, likeness_table = "BLOSSUM62"):
        """
        Instantiation of the aligner object. It loads needed data.
        
        @param cc: A coiled coil object containing score and heptad data. If none is 
                given, then it defaults to a standard one.
        @type cc: CoiledCoil
        @param likeness_table: Name of the file used for loading the residues likeness table.
        @type likeness_table: string
        """
        self.like_scores, self.charge_scores = self.createTables(likeness_table)
        self.mycc = cc or CoiledCoil()
        self.chain_alignments = {}
        self.reg_alignments = {}
        
    
    def createTables(self,table = ""):
        """
        Load residue likelyhood scores for BLOSSUM 62 table.
        Creates from scratch a table for residue charge comparison.
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
        
        self.chain_alignments["heptads"] = alignf(a_s,b_s,scoreAlignFun)
        scc = max(self.chain_alignments["heptads"])
        
        ## Chain similarity
        self.chain_alignments["res_like"] = alignd(a, b, self.like_scores)
        chs = max(self.chain_alignments["res_like"])
        
        ## Charges similarity
        self.chain_alignments["charges"] = alignd(a, b, self.charge_scores)
        crs = max(self.chain_alignments["charges"])
        
        return {"heptads": scc , "res_like": chs , "charges":  crs}
    

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
        scores_a = flatten(a,"abcdefg"*total_regs_a,self.mycc)
        scores_b = flatten(b,"abcdefg"*total_regs_b,self.mycc)
        
        a_s,b_s = scores2String(scores_a,scores_b)
        
        scc = alignf(a_s,b_s,scoreAlignFun)
        
        new_scc = []
        for t in scc:
            new_scc.append(( t[0], start_a+(t[1]*7) ))
        
        scc = ( max(scc)[0], start_a+max(scc)[1]*7)
        
        self.reg_alignments["heptads"]= new_scc
        
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
        
        self.reg_alignments["charges"]=subvals_charge
        self.reg_alignments["res_like"]=subvals_res
        
        crs = max(subvals_charge)
        chs = max(subvals_res)
       
        #~ return {"heptad": scc, "charges": max(subvals_charge),"res_like": max(subvals_res)}
        a = ch_a;b=ch_b;reg_a=bra;reg_b=brb
        return {"heptads":scc , "res_like":chs, "charges": crs } ##(self.getFitnessScore("charges",a,b,crs[1],reg_a,reg_b),crs[1])}
    
    def normalizeScores(self):
        """
        Normalizes all the alignment scores to a 0-1 range.
        """
        for i in self.chain_alignments:
            mymax = self.chain_alignments[i][0]
            mymin = self.chain_alignments[i][0]
            if len(self.chain_alignments[i]) >1:
                for j in self.chain_alignments[i]:
                    mymax = max(mymax,j)
                    mymin = min(mymin,j)
               
                for k in range(len(self.chain_alignments[i])):
                    self.chain_alignments[i][k] = ( (self.chain_alignments[i][k][0]-mymin[0]) / (mymax[0]-mymin[0]) , self.chain_alignments[i][k][1])
        
        for i in self.reg_alignments:
            mymax = self.reg_alignments[i][0]
            mymin = self.reg_alignments[i][0]
            if len(self.reg_alignments[i]) >1:
                for j in self.reg_alignments[i]:
                    mymax = max(mymax,j)
                    mymin = min(mymin,j)
               
                for k in range(len(self.reg_alignments[i])):
                    self.reg_alignments[i][k] = ( (self.reg_alignments[i][k][0]-mymin[0]) / (mymax[0]-mymin[0]) , self.reg_alignments[i][k][1])

    def copy(self):
        """
        Returns a copy where tables are referenced (and then not loaded again).
        For speeding up purposes.
        """
        ca = CoiledAlign(self.mycc)
        ca.like_scores = self.like_scores
        ca.charge_scores = self.charge_scores
        
        return ca
        
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
        self.assertEqual(res,{'res_like': (18.0, 7), 'charges': (2.6000000000000001, 7), 'heptads': (404, 7)})
    
    
    def test_align_regs(self):
        """Best alignment for registers test cases"""
        
        res = self.ca.alignRegisters("RRRLLLLLLLRRRRRRLLLLLLLRRR", "LLLRRRR","efgabcdefgabcdefgabcdefgab","abcdefg")
        self.assertEqual(res["charges"],(2, 10))
        res = self.ca.alignRegisters("RRRLLLLLLLRRRRRRLLLLLLLRRR", "RRRLLLL","efgabcdefgabcdefgabcdefgab","abcdefg")
        self.assertEqual(res["charges"],(3.20, 10))
        res = self.ca.alignRegisters("RRRLLLLLLLRRRRRRLLLLLLLRRR", "RRRRRRR","efgabcdefgabcdefgabcdefgab","abcdefg")
        self.assertEqual(res["charges"],(11, 10))
    
    def test_normalize(self):
        """ Data normalization test_case """
        data = {}
        data["heptad"] =[ (1.,0),(3.,0),(10.,0), (-1.,0), (5.,0)]
        data["charges"] =[ (1.,0),(3.,0),(10.,0), (-1.,0), (5.,0)]
        data["res_like"] =[ (1.,0),(3.,0),(10.,0), (-1.,0), (5.,0)]
        
        data2 = {}
        data2["heptad"] =[ (1.,0)]
        data2["charges"] =[ (1.,0)]
        data2["res_like"] =[ (1.,0)]
        
        self.ca.chain_alignments = data
        self.ca.reg_alignments = data2
        
        self.ca.normalizeScores()
        self.assertEqual( self.ca.chain_alignments["heptad"][0][0], (1.-(-1)) / (10-(-1))) 
    
    def test_copy(self):
        """Copy test case"""
        ca2 = CoiledAlign(self.cc)
        ca3 = self.ca.copy()
        
        self.assertEqual( ca3.like_scores[("B","N")],3)
    
if __name__ == '__main__':
    BT.localTest()    
    