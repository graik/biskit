
from Biskit.PDBModel import PDBModel
from numpy import compress,transpose,where,cross,array,greater,matrix,array,greater
import Biskit.molUtils as MU 
from emath import norm,sphericalAngles
import cPickle 
import os
from resTools import doReorientation

class FoldDataCreator:
    
    ##TODO: Usar IA o simples comparaciones de una BD
    
    def __init__(self,db="./fold_db"):
        self.dbpath = db
        self.res_char = {}
        self.inner_radius = 5
        self.outer_radius = 10
        self.segment = 1
        self.inc_angle = 45
        self.constant_vol = False
        self.wantInfoRate = False
    
    def process_list (self,list):
        pass
    
    def save (self, where):
        """
        Saves the current data stored in the class to disk.
        """
        try:
            f = open (where,"wb")
            
            cPickle.dump(self.res_char,f,cPickle.HIGHEST_PROTOCOL)

        except (cPickle.UnpicklingError, IOError,EOFError):
            print "ERROR!!!"
        
        f.close()
        
        
    def restore (self, where):
        """
        Loads one or multiple files from disk.
        """
        try:
            f = open (where,"rb")
            
            self.res_char= cPickle.load(f)

        except (cPickle.UnpicklingError, IOError,EOFError):
            print "ERROR!!!"
        
        f.close()
        
        
    def process_protein (self, mymodel = None, mychain = 0):
        if not isinstance(mymodel,PDBModel):
            model = PDBModel(mymodel)
        
        raw = model.clone()
        
        ## Clean the structure a bit
        model = model.compress( model.maskProtein())

        model.report()
        
        try:
            chain = model.takeChains([mychain])
            raw = model.takeChains([mychain])
        except:
            return
        
        ## We are only going to work with N atoms
        chain = chain.compress(chain.maskFrom( 'name', 'N' ))
        
        for i in range(1,len(chain.xyz)-1) :
            
            key = (MU.singleAA([chain.atoms["residue_name"][i-1]])[0],MU.singleAA( [chain.atoms["residue_name"][i]])[0],MU.singleAA([chain.atoms["residue_name"][i+1]])[0])
            print i,": ",chain.atoms["residue_name"][i]
            
            env = self.getEnvironment( raw, key, i, self.inner_radius,self.outer_radius,self.segment,self.inc_angle,self.constant_vol ) 
            
            if not key in self.res_char:
                self.res_char[key] = []
            
            self.res_char[key].append(env)
                
            
    def getEnvironment (self, mymodel=None,env_res = [], res = 0, inradius = 5,outradius=10,segment = 1, incangle= 45 ,constanvol = False):
        ## TODO: Offer 2 options, constant volume or not (default)
        ## save also probabilities for each component in a box 
        
        ## Reorient the whole molecule
        model = mymodel.clone()
        
        model = doReorientation(mymodel)
        
        ## Get the sphere data for a sphere of radius d
        ## All interesting points inside the sphere centered in Ca 
        residue = model.resModels()[res]
        center = residue.xyz[1]
        
        selected = []
        for i in range(len(model.xyz)):
            def_vect = (model.xyz[i] - center)
            n = norm(def_vect)
            if  n >= inradius and n <= outradius:
                a,b = sphericalAngles(def_vect)
                selected.append((model['residue_name'][i],a*180/3.14159,b*180/3.14159,n))
        
        len_sel = len(selected)
        
        
        ## Calculate information by res
        probs = {}
        for s in selected:
            if s[0] in probs:
                probs[s[0]] += 1 
            else:
                probs[s[0]] = 1
        
        for k in probs.keys():
            if k in env_res:
                probs[k] = 1./len_sel
            else:
                probs[k] = 1./probs[k]
       
        
        
        ## Then get the environment vector
        penv = [([([{}]* (360/incangle))]*(360/incangle))]*((outradius-inradius)/segment)
        env = [([([["WAT"]]* (360/incangle))]*(360/incangle))]*((outradius-inradius)/segment)
        
        len_env = (360/incangle)**2 *(outradius-inradius)/segment
        

        for s in selected:
            a = int(s[3]-inradius)
            b = int(s[1]/incangle)
            c = int(s[2]/incangle)
               
            if not s[0] in penv[a][b][c]:
                penv[a][b][c][s[0]] = probs[s[0]]
            else:
                penv[a][b][c][s[0]] += probs[s[0]]
        
            env[a][b][c]= self.best(penv[a][b][c])
        
        if self.wantInfoRate == True:
            return penv
        else:
            return env
        
    def best(self,d={}):
        ## retorna el mejor de un diccionario de puntuacion de info
        totals = []
        for k in d.keys():
            totals.append((d[k],k))
        
        return max(totals)[1]
        
##############
## Test
##############
import Biskit.test as BT
import Biskit.tools as T


class Test(BT.BiskitTest):
    """ Test cases for fold data creation"""

    def prepare(self):
        self.f = FoldDataCreator()
        

    def cleanUp( self ):
        
        pass
    
    def test_DB_Creation(self):
        """DB Creation test cases"""
        self.f.process_protein(T.testRoot()+"/polysys/FAKE.pdb")
        #~ self.f.process_protein(T.testRoot()+"/polysys/FAKE.pdb")
    
    def test_Best(self):
        """ test Best function """
        d = {"LOL":4.0,"LEL":4.2,"LIL":3.8}
        print self.f.best(d)
        
    def test_Save_and_Restore(self):
        """ test Save and Restore functions """
        self.f.process_protein(T.testRoot()+"/polysys/FAKE.pdb")
        self.f.save("./this")
        f2 = FoldDataCreator()
        f2.restore("./this")
        self.assertEqual(self.f.res_char,f2.res_char)
        os.remove("./this")
        
        
if __name__ == '__main__':
    BT.localTest()    
