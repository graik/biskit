
from Biskit.PDBModel import PDBModel
from numpy import compress,transpose,where,cross,array,greater,matrix,array,greater
import Biskit.molUtils as MU 
from emath import norm,sphericalAngles

class FoldDataCreator:
    
    ##TODO: Usar IA o simples comparaciones de una BD
    
    def __init__(self,db="./fold_db"):
        self.dbpath = db
        self.res_char = {}
    
    def process_list (self,list):
        pass
    
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
            print i
            env = self.getEnvironment( raw, i ) 
            
            if not key in self.res_char:
                self.res_char[key] = []
            
            self.res_char[key] = []
                
            
    def getEnvironment (self, model=None, res = 0, inradius = 5,outradius=10,segment = 1, incangle= 45 ,constanvol = False):
        ## TODO: Offer 2 options, constant volume or not (default)
        ## save also probabilities for each component in a box 
        
        ## Reorient the whole molecule
        pass
        
        ## Get the sphere data for a sphere of radius d
        
        ## All interesting points inside the sphere centered in Ca 
        residue = model.resModels()[res]
        center = residue.xyz[1]
        this =  residue.atoms['residue_name'][0]
        
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
            probs[k] = 1./probs[k]
        
        
        ## Then get the environment vector
        penv = [([([{}]* (360/incangle))]*(360/incangle))]*((outradius-inradius)/segment)
        env = [([([["WAT"]]* (360/incangle))]*(360/incangle))]*((outradius-inradius)/segment)
        
        len_env = (360/incangle)**2 *(outradius-inradius)/segment
        if this in probs:
            probs[this] = 1 / len_sel
            
         
        for s in selected:
            
               
            if not s[0] in penv[int(s[3]-inradius)][int(s[1]/incangle)][int(s[2]/incangle)]:
                penv[int(s[3]-inradius)][int(s[1]/incangle)][int(s[2]/incangle)][s[0]] = probs[s[0]]
                
               
            else:
                penv[int(s[3]-inradius)][int(s[1]/incangle)][int(s[2]/incangle)][s[0]] += probs[s[0]]
        
            env[int(s[3]-inradius)][int(s[1]/incangle)][int(s[2]/incangle)]= best(penv[int(s[3]-inradius)][int(s[1]/incangle)][int(s[2]/incangle)][s[0]])
        
        return env
        
    def best(self,d={}):
        ## retorna el mejor de un diccionario de puntuacion de info
        pass
        
##############
## Test
##############
import Biskit.test as BT
import Biskit.tools as T


class Test(BT.BiskitTest):
    """ Test cases for Polyfret"""

    def prepare(self):
        pass

    def cleanUp( self ):
        #~ os.system("rm -rf "+T.testRoot()+"/polysys/residues_db")
        pass
    
    def test_DB_Creation(self):
        """DB Creation test cases"""
        f = FoldDataCreator()
        f.process_protein(T.testRoot()+"/polysys/FAKE.pdb")
        

if __name__ == '__main__':
    BT.localTest()    
