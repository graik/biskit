
from Biskit.PDBModel import PDBModel
from numpy import compress,transpose,where,cross,array,greater,matrix,array,greater
import Biskit.molUtils as MU 
from vectors import norm,sphericalAngles
import cPickle 
import os
from restools import doAAReorientation,aalist2intlist
import Biskit.tools as T
from pyfann import libfann
import random

class FoldDataCreator:
    
    ##TODO: Usar IA o simples comparaciones de una BD
    
    def __init__(self,db="./fold_db"):
        self.dbpath = db
        ## Try to create the folder
        if( not os.access(self.dbpath, os.F_OK)):
            os.mkdir(self.dbpath)
            
        self.res_char = {}
        self.res_trans = {}
        self.inner_radius = 5
        self.outer_radius = 15
        self.segment = 2
        self.inc_angle = 45
        self.constant_vol = False
        self.wantInfoRate = False
        self.len_env = (360/self.inc_angle)**2 *(self.outer_radius-self.inner_radius)/self.segment
        
        ## Try to load a trained Neural Network or if not, create a new one otherwise
        self.nn = libfann.neural_net()
        
        
        self.restore(self.dbpath)
        
      
    def save (self, where=""):
        """
        Saves the current data stored in the class to disk.
        """
        
        try:
            if where == "":
                where = self.dbpath
            
            if( not os.access(where, os.F_OK)):
                os.mkdir(where)
                
            f = open (where+"/res_env","wb")
            print "-Saving residue environmental info to",where+"/res_env"
            cPickle.dump(self.res_char,f,cPickle.HIGHEST_PROTOCOL)
            
            f.close()
            
            self.nn.save(self.dbpath+'/neural.network')
            print "-Saving neuronal network to",where+"/neural.network"
            
            
        except (cPickle.UnpicklingError, IOError,EOFError):
            print "ERROR WHILE SAVING"
        

    def restore (self, where=""):
        """
        Loads one or multiple files from disk.
        """
        if where == "":
                where = self.dbpath
        
        print "- Trying to open neural network in:", self.dbpath+"/neural.network"
        
        if not self.nn.create_from_file(self.dbpath+"/neural.network"):
            print "- Impossible.Creating a new one."
            ## Init fann neural network (based on example)
            connection_rate = 0.8
            learning_rate = 0.5 ## This is the default value
            num_input = self.len_env + 3
            num_neurons_hidden = num_input
            num_output = 3

            desired_error = 0.0001
            max_iterations = 1000000
            iterations_between_reports = 1000

            self.nn.create_sparse_array(connection_rate,(num_input, num_neurons_hidden/2,num_neurons_hidden/4,num_neurons_hidden/8 ,num_output))
            self.nn.randomize_weights(-0.77,0.77)
            self.nn.set_learning_rate(learning_rate)
            self.nn.set_training_algorithm(libfann.TRAIN_RPROP)
            self.nn.set_train_error_function(libfann.ERRORFUNC_TANH)
            self.nn.set_activation_function_output(libfann.LINEAR)
            self.nn.set_activation_function_layer(libfann.SIGMOID,1)
            self.nn.set_activation_function_layer(libfann.SIGMOID_SYMMETRIC,2)
            self.nn.set_activation_function_layer(libfann.SIGMOID_SYMMETRIC,3)
            
            #~ self.nn.set_activation_function_hidden(libfann.SIGMOID_SYMMETRIC)
            print "Total neurons in network:",self.nn.get_total_neurons()
        else:
            print "OK"
        
        try:
            
            print "-Trying to load residue environmental info from:",where+"/res_env"
            
            if os.access(where+"/res_env", os.F_OK):
                f = open (where+"/res_env","rb")
                
                self.res_char= cPickle.load(f)
                
                f.close()
                print "OK"
        except (cPickle.UnpicklingError, IOError,EOFError):
            print "ERROR WHILE RESTORING"
        
        
        
        
    def process_protein (self, mymodel = None, mychain = 0):
        if not isinstance(mymodel,PDBModel):
            model = PDBModel(mymodel)
        else:
            model = mymodel
        raw = model.clone()
        
        ## Clean the structure a bit
        model = model.compress( model.maskProtein())
        
        
        try:
            chain = model.takeChains([mychain])
            raw = chain.clone()
        except:
            return
            
        residuos = chain.resModels()
        ## We are only going to work with N atoms
        chain = chain.compress(chain.maskFrom( 'name', 'N' ))
        
        for i in range(1,len(residuos)-1) :
            
            key = (residuos[i-1].sequence(),residuos[i].sequence(),residuos[i+1].sequence())
            print i,": ",residuos[i].sequence()
            
            env = self.getEnvironmentalData( residuos[i],raw, key, i, self.inner_radius,self.outer_radius,self.segment,self.inc_angle,self.constant_vol ) 
             
            if not key in self.res_char:
                self.res_char[key] = []
            #~ else:
                #~ print "repeated:",key ## Depending on a comparison put the regor not
            
            
            self.res_char[key].append(env)
        
        del residuos 
        del model
        del raw
            
    def getEnvironmentalData (self, resm=None,mymodel=None,env_res = [], res = 0, inradius = 5,outradius=10,segment = 1, incangle= 45 ,constanvol = False):
        ## TODO: Offer 2 options, constant volume or not (default)
        ## save also probabilities for each component in a box 
        
        
        ## Reorient the residue
        model = mymodel.clone()
        recenter = [ resm.xyz[0][0],resm.xyz[0][1],resm.xyz[0][2]]## Make  this propperly
        residue, R = doAAReorientation(resm.clone())
        
        ## center model
        model.xyz = model.xyz - recenter 
        del recenter
        
        ## rotate model
        for i in range(len(model.xyz)):
            model.xyz[i] = model.xyz[i] * R 
        
        residues = model.resModels()
        ## Translation vector is the vector for going to the next N
        ## so if we are centered is just the next one  ^.^
        trans = [residues[res+1].xyz[0][0],residues[res+1].xyz[0][1],residues[res+1].xyz[0][2]]
        del residues
        
        ## Get the sphere data for a sphere of radius d
        ## All interesting points inside the sphere centered in Ca 
        center = residue.xyz[1] ## Alpha Carbon TODO choose atom propperly
        
        selected = []
        for i in range(len(model.xyz)):
            def_vect = (model.xyz[i] - center)
            n = norm(def_vect)
            if  n >= inradius and n <= outradius:
                a,b = sphericalAngles(def_vect)
                if a<0.:
                    a = 2*3.1415 + a
                if a>=2*3.141:
                    a = 2*3.141
                if b<0.:
                    b = 2*3.1415 + b
                if b>=2*3.141:
                    b = 2*3.141
                selected.append((model['residue_name'][i],a*180/3.14159,b*180/3.14159,n))
        
        len_sel = len(selected)
        del model
        
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
        penv = []
        env  = []
        for i in range((outradius-inradius)/segment):
            env.append([])
            penv.append([])
            for j in range(360/incangle):
                env[i].append([])
                penv[i].append([])
                for k in range(360/incangle):
                    env[i][j].append('WAT')
                    penv[i][j].append({})
                    
               
        len_env = (360/incangle)**2 *(outradius-inradius)/segment
        
        
        for s in selected:
            a = int((s[3]-inradius)/segment)
            b = int(s[1]/incangle)
            c = int(s[2]/incangle)
            
            if not s[0] in penv[a][b][c]:
                penv[a][b][c][s[0]] = probs[s[0]]
            else:
                penv[a][b][c][s[0]] += probs[s[0]]
        
            env[a][b][c]= self.best(penv[a][b][c])
        
        del probs
        del selected
        
        if self.wantInfoRate == True:
            del env
            return penv, trans
        else:
            del penv
            return env , trans
        
    def best(self,d={}):
        ## Retorna el mejor de un diccionario de puntuacion de info
        totals = []
        for k in d.keys():
            totals.append((d[k],k))
        
        return max(totals)[1]
        
    def consumeData(self):
        aas = MU.allAA()
        cor ={}
        for i in range(len(aas)):
            cor[aas[i]]=i
        del aas
        
        environments=[]
        outputs = []

        for key in self.res_char.keys():
            for l in self.res_char[key]:
                env = []
                for i in range(len(l[0])):
                    for j in range(len(l[0][i])):
                        env = env + aalist2intlist(l[0][i][j])
                
                env = [cor[key[0]]+25,cor[key[1]]+25,cor[key[2]]+25]+env
                environments.append(env)
                outputs.append(l[1])
                
        
        self.genDataFile("tempdata",environments,outputs)
        trainingdata= libfann.training_data()
        trainingdata.read_train_from_file("tempdata")
        
        self.nn.train_on_data(trainingdata,500,1,0.001)
        
         
        del cor
        del self.res_char
        del environments
        del outputs
        self.res_char = {}
        
    def genDataFile(self,path,input,output):
        file = open(path,"w")
        
        file.write("%d %d %d"%(len(input),len(input[0]),len(output[0])))
        
        
        for i in range(len(input)):
            file.write("\n")
            for j in range(len(input[i])):
                file.write("%f "%(input[i][j]))
            file.write("\n")
            for j in range(len(output[i])):
                file.write("%f "%(output[i][j]))
            
        file.close()
    
    def getTranslations(self,mymodel = None, mychain = 0):
        aas = MU.allAA()
        cor ={}
        for i in range(len(aas)):
            cor[aas[i]]=i
        del aas
        
        
        if not isinstance(mymodel,PDBModel):
            model = PDBModel(mymodel)
        else:
            model = mymodel
            
        points = []
        
        ## Clean the structure a bit
        model = model.compress( model.maskProtein())

        try:
            chain = model.takeChains([mychain])
            raw = chain.clone()
        except:
            print "ERROR while extracting chain"
            return
        print "Generating..."
        residuos = chain.resModels()
        
        
        ## We are only going to work with N atoms
        chain = chain.compress(chain.maskFrom( 'name', 'N' ))
        
        ## Random translation for the first
        points.append(array([(random.random()*2)-1+2.65,(random.random()*2)-1+2.65,(random.random()*2)-1+2.65]))
        
        print "Hay tantos residuos",len(residuos)
        for i in range(1,len(residuos)-1) :
            
            key = (residuos[i-1].sequence(),residuos[i].sequence(),residuos[i+1].sequence())
            print i,": ",residuos[i].sequence()
            
            environment = self.getEnvironmentalData( residuos[i],raw, key, i, self.inner_radius,self.outer_radius,self.segment,self.inc_angle,self.constant_vol ) 
            
            env = []
            for i in range(len(environment[0])):
                for j in range(len(environment[0][i])):
                    env = env + aalist2intlist(environment[0][i][j])
            
            env = [cor[key[0]]+25,cor[key[1]]+25,cor[key[2]]+25]+env
            
            points.append(array(self.nn.run(env)))
        
        del residuos 
        del model
        del raw
        
        
        ## add random translation for the last residue
        points.append(array([(random.random()*2)-1+2.65,(random.random()*2)-1+2.65,(random.random()*2)-1+2.65]))
        
        
        return array(points)
##############
## Test
##############
import Biskit.test as BT



class Test(BT.BiskitTest):
    """ Test cases for fold data creation"""

    def prepare(self):
        self.f = FoldDataCreator(T.testRoot()+"/polysys/this")
        

    def cleanUp( self ):
        pass
    
    def test_DB_Creation(self):
        """DB Creation test cases"""
        self.f.process_protein(T.testRoot()+"/polysys/FAKE.pdb")
        
        ## Those two are the same protein in two orientation.May have the same
        ## environment vectors
        self.f.process_protein(T.testRoot()+"/polysys/envtesta.pdb")
        self.f.process_protein(T.testRoot()+"/polysys/envtestb.pdb")
        
        ## Train the network with the first data
        self.f.consumeData()
        
        ## Process a big protein
        self.f.process_protein(T.testRoot()+"/polysys/1HUY.pdb")
        ## And train!
        self.f.consumeData()
        
    def test_Best(self):
        """ test Best function """
        d = {"LOL":4.0,"LEL":4.2,"LIL":3.8}
        self.assertEqual( self.f.best(d),"LEL")
        
    def test_Save_and_Restore(self):
        """ test Save and Restore functions """
        print
        self.f.process_protein(T.testRoot()+"/polysys/FAKE.pdb")
        print 
        print
        self.f.save(T.testRoot()+"/polysys/this")
        print 
        print 
        f2 = FoldDataCreator(T.testRoot()+"/polysys/this")
        self.assertEqual(self.f.res_char,f2.res_char)
        f2.consumeData()
        
        
if __name__ == '__main__':
    BT.localTest()    
