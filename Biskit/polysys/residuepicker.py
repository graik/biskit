from Biskit.PDBModel import PDBModel
from numpy import compress,transpose,where,cross,array,matrix,array,cross, dot
import cPickle 
import os
from math import acos,sqrt,atan2,pi
import emath 

from quaternion import rotquat, rotmat
import Biskit.molUtils as MU 
from resTools import doAAReorientation, calcPlane, genPoints,planarize
from tools import lendepth
import random

from vectors import angle,add, sub

class residuePicker:
    
    def __init__(self,databasepath = "./residues_db/",verbose = False):
        """
        Initializes dictionaries and lists used by the residue picker:
        
        - index - Is a dictionary of lists which stores the relative path of stored residues without extension: e.g. "D/D0" for
                    the first Asp residue stored.
        
        - apvect - Also a dictionary for storing usefull data about the index residues.
        
        - last_res - Dictionary of the last number used for indexing db residues (for naming purposes).
        
        - reslib - Dictionary of lists with all the PDBModels for each residue in db.

        - clusters - Dictionary of lists of lists tracking backbone similarities.
        
        @param databasepath: Path of the db folder that is going to be used.
        @type databasepath: string
        @param verbose: If True the screen fills with lots of usefull text.
        @type verbose: bool
        """
        self.verbose = verbose
        self.dbpath = databasepath
        
        try:
            f = open (self.dbpath+'index.db',"rb")
            self.index = cPickle.load(f)

        except (cPickle.UnpicklingError, IOError,EOFError):
  
            if( not os.access(self.dbpath, os.F_OK)):
                os.mkdir(self.dbpath)
            f = open (self.dbpath+'index.db',"wb")
            self.index = {}
        f.close()
        
        try:
            f = open (self.dbpath+'apvect.db',"rb")
            self.apvect = cPickle.load(f)

        except (cPickle.UnpicklingError, IOError,EOFError):
            f = open (self.dbpath+'apvect.db',"wb")
            self.apvect = {}
        f.close()
        
        try:
            f = open (self.dbpath+'clusters.db',"rb")
            self.clusters = cPickle.load(f)

        except (cPickle.UnpicklingError, IOError,EOFError):
            f = open (self.dbpath+'clusters.db',"wb")
            self.clusters = {}
        f.close()
        
        ## Create directory structure
        for k in MU.allAACodes():
            if not os.access(databasepath+k, os.F_OK):
                os.mkdir(databasepath+k)
        
        self.last_res = {}
        self.reslib = {}
        ## Initialize them
        for k in MU.allAACodes():
            if not k in self.apvect.keys():
                self.apvect[k] = []
            if not k in self.index.keys():
                self.index[k] = []
            self.last_res[k] = 0
            
    
        ## Then see what the last used id was for each residue type
        
        for k in self.index.keys():
            if self.index[k]!=[]:
                self.last_res[k] = int(self.index[k][-1][1:])
    
        ## Create reslib
       
        ## load all residues inside it
        for k in self.index.keys():
            for i in self.index[k]:
                self.reslib[i]=PDBModel(self.dbpath+i[0]+"/"+i+".pdb")
        
        ## create the cluster
        for k in MU.allAACodes():
            if not k in self.clusters:
                self.clusters[k]=[]
    
    
        
    def save(self):
        """
        Saves the index, clusters and derived data of the dabase.
        """
        
        f = open (self.dbpath+'index.db',"wb")
        cPickle.dump(self.index,f,cPickle.HIGHEST_PROTOCOL)
        f.close()
        
        if self.verbose:
            for k in self.index:
                print k, len(self.index[k])
        
        f = open (self.dbpath+'apvect.db',"wb")
        cPickle.dump(self.apvect,f,cPickle.HIGHEST_PROTOCOL)
        f.close()
        
        f = open (self.dbpath+'clusters.db',"wb")
        cPickle.dump(self.clusters,f,cPickle.HIGHEST_PROTOCOL)
        f.close()
        if self.verbose:
            print "Clusters: ",self.clusters

    def extractFromPDB ( self,path = "" ):
        """
        
        """
        ## Clean the structure a bit
        model = PDBModel(path)
        model = model.compress( model.maskProtein())
        
        ## Also remove all hidrogens
        model = model.compress(model.maskHeavy())
        
        if self.verbose:
            model.report()
            
        for c in range( model.lenChains() ):
            
            chain = model.takeChains( [c] )
            
            residuos = chain.resModels()
            
            for i in range(0,len(residuos)-1) :
                if self.verbose:
                    print "-"+ residuos[i].sequence() + str(self.last_res[residuos[i].sequence()])
                r = residuos[i].sequence()
                residue = residuos[i]
                myindex = r+str(self.last_res[r])
                
                ## Fetch neighbouring residues
                nextres = residuos[i+1]
                rnext = nextres.sequence()
                
                    
                ## First purge: it has all the needed atoms?
                if len(residue.xyz) == len(MU.aaAtoms[MU.single2longAA(r)[0]])-1:
                    
                    ## Hunt for N & C to reorient the aa
                    int_atoms = ['N','CA','C']
                    sel_atoms= where(residue.maskFrom( 'name', int_atoms ))
                    sel_atoms_next= where(nextres.maskFrom( 'name', int_atoms ))
                    Ni  = residue.xyz[sel_atoms[0][0]]
                    Ca = residue.xyz[sel_atoms[0][1]]
                    C  = residue.xyz[sel_atoms[0][2]]
                    Ninext  = nextres.xyz[sel_atoms_next[0][0]]
                    Canext = nextres.xyz[sel_atoms_next[0][1]]
                    Cnext  = nextres.xyz[sel_atoms_next[0][2]]
    
                    
                    Translation = Ninext-Ni
                    Translation = [Translation[0],Translation[1],Translation[2]]
                    Plane = [(C-Ni),(Ca-Ni)]
                    Plane_next=[(Cnext-Ninext),(Canext-Ninext)]
                    residue , R , rv = doAAReorientation(residue)
                    
                    ## Convert the obtained vectors
                    Translation = Translation*R
                    Translation = [Translation[0,0],Translation[0,1],Translation[0,2]]
                    Plane[0]= Plane[0]*R
                    Plane[0] = [Plane[0][0,0],Plane[0][0,1],Plane[0][0,2]]
                    Plane[1]= Plane[1]*R
                    Plane[1] = [Plane[1][0,0],Plane[1][0,1],Plane[1][0,2]]
                    Plane_next[0]= Plane_next[0]*R
                    Plane_next[0] = [Plane_next[0][0,0],Plane_next[0][0,1],Plane_next[0][0,2]]
                    Plane_next[1]= Plane_next[1]*R
                    Plane_next[1] = [Plane_next[1][0,0],Plane_next[1][0,1],Plane_next[1][0,2]]
                    
                    Plane_angle = angle(cross(Plane[0],Plane[1]),cross(Plane_next[0],Plane_next[1]))
                    add_it = False
                    foundback = False
                    exit = False
                    
                    ## Second purge: RMSD and clusterization
                    for group in self.clusters[r]:
                            ## We compare with the first
                            (bb,sd) = self.areTheSame(residue,self.reslib[group[0]])
                            
                            ## If it has the same backbone than the group head, then we have to 
                            ## find further similarities to se if it is already in the list.
                            foundsd = False
                           
                            
                            if bb and not exit:
                                
                                foundback = True
                                
                                mycluster = group[0]
                                
                                exit = True ## if we have found the bck then the next step is to search for sd inside
                                            ## this group
                                
                                for m in group:
                                    (bb,sd) = self.areTheSame(residue,self.reslib[m])
                                    if sd:
                                       foundsd= True 
                                       similarto = m 
                                
                                if foundsd:
                                    ## If backbone and sidechains are identical we don't need to 
                                    ## add the residue
                                    addit=False
                                    if self.verbose:
                                        print "Too similar to other ("+ similarto +")"
                                    myindex = similarto
                                    
                                else:##if not foundsd:
                                    group.append(myindex)
                                    add_it = True
                                    if self.verbose:
                                        print "Backbone found but different sidechain:adding"
                                    
                    if not foundback:
                        ## We haven't found the backbone so we need to create a new
                        ## group.
                        self.clusters[r].append([myindex])
                        mycluster = myindex
                        add_it = True
                        if self.verbose:
                            print "New group (backbone not found)"
                    

                    if add_it:
                        ## Save the residue
                        residue.renumberResidues()
                        residue['serial_number'] = range(1,residue.lenAtoms()+1)
                        residue.writePdb(self.dbpath+r+"/"+myindex+".pdb")
                        
                        ## Add and index it
                        self.reslib[myindex]=residue
                        self.index[r].append(myindex)
                        self.last_res[r]+=1
                    
                    
                    key = (r,rnext)
                    new_entry =  {'index':myindex,'next':rnext,'cluster':mycluster,'translation':Translation,'plane_next':Plane_next,'planes_angle':Plane_angle}
                    
                    if not key in self.apvect.keys():
                        ## if the key isn't we just add the new register
                        self.apvect[key] = [ new_entry ]   
                    else:
                        ## else we have to see if there's something similar inside
                        foundsimilar = False
                        for i in self.apvect[key]:
                            #~ print new_entry,i
                            if angle(new_entry['translation'],i['translation']) < 0.06 and abs(new_entry['planes_angle']-i['planes_angle']) < 0.06:
                                foundsimilar = True
                                print "Found Similar"
                        if foundsimilar == False:
                            self.apvect[key].append(new_entry)
                else:
                    if self.verbose:
                        print "Diferent number of atoms."
        self.save()
    
    
        
    def areTheSame(self,a,b,bckbthres = 0.2, sdthres= 0.3):
        """
        It compares two residues to see if they are too similar to be stored in the database.
        First it compares backbone, and then it compares the sidechain. One can choose the 
        thresholds used for each comparison, but default ones are good enough.
        
        @param a: Residue to be compared.
        @type a: PDBModel
        @param b: Residue to be compared. 
        @type b: PDBModel
        @param bckbthres: Threshold for backbone comparisons.
        @type bckbthres: float
        @param sdthres: Threshold for sidechain comparisons.
        @type sdthres: float
        
        @return: True if the too residues are very similar in conformation.
        @rtype: bool
        """
        
        bcbsim = False
        sdsim = False
        
       
        abk_mask = a.maskFrom( 'name', ["N","CA","C","O"])
        bbk_mask = b.maskFrom( 'name', ["N","CA","C","O"])
        
        abk = a.take(where(abk_mask)[0])
        bbk = b.take(where(bbk_mask)[0])
        
        
        
        if abk.rms(bbk)<= bckbthres:
           bcbsim = True
           
        asd = a.take(range(4,len(b['element'])))
        bsd = b.take(range(4,len(b['element'])))
        if a.sequence()=='G':
            sdsim = True
        else:
            if  asd.rms(bsd)<= sdthres:
                sdsim = True
                
        
        return bcbsim , sdsim

    
    def randomRes(self, res = 'A'):
        """
        Returns a random residue loaded from the library of type resname.
        
        @param resname: Single letter code of aa.
        @type resname: string
        
        @return: A residue model from the database.
        @rtype: PDBModel
        """
        if isinstance(res,list):
            choosen = res[random.randrange(len(res))]
        else:
            choosen = self.index[res][random.randrange(len(self.index[res]))]
       
    
        return choosen
    
    def chooseRes(self,r='A',rnext = 'A',reg = {},simfun=None):
        candidates = self.apvect[(r,rnext)]
        
        sim = []
        
        for i in range(len(candidates)): # hacerlo por indices
            sim.append((simfun(reg,candidates[i]),i))
        
         
        #####################################
        ##### TODOOOOOOOOOOOOOOOOOOOOOOOOO ##
        #####################################
        ## devolver  randomres de su cluster? Dar la opcion (o elegir segun steric clashes)
            
        return candidates[max(sim)[1]]
    
    def defaultSimFun(self,reg1,reg2):
        return angle(reg1['translation'],reg2['translation']) / pi
        
    
    
    def createChain(self,seq="AA",points=[[0,0,0]],simfun = None):
        print "\n\n\n\n\n"
        T=array([0,0,0])
        if simfun == None:
            simfun = self.defaultSimFun
        
        r = seq[0]
        chain = PDBModel(self.dbpath+r+"/"+self.randomRes(r)+".pdb")
        
        ## Suponiendo siempre que en xyx[0] tenemos N 
        ## TODO no suponerlo
        N = chain.xyz[0]
        i = 0
        print "points ",len(points)
        for c in seq[1:len(seq)]:
            rnext = c
            print i
            reg ={"translation":list(N - points[i])}
            nextres_reg = self.chooseRes(r,rnext,reg,simfun)
            print "Key: (",r,",",rnext,") Next:", nextres_reg
            
            ## Cargarlo como PDBModel
            ## !!!!!!!!!!!!!!!! Nohay manera de saber el ID completo de next
            ## a no ser que se haga el procesado en dos tiradas.La segunda se conoce
            ## las IDs propias o del representante
            res = PDBModel(self.dbpath+rnext+"/"+self.randomRes(nextres_reg['next'])+".pdb")
            
            ## Make rotationmatrix from this plane to the other
            ## Generar una matriz que rote el otro plano al plano inicial (xz)
            ## Y multiplicarlo todo por esa matriz
            Nlast = chain.compress(chain.maskFrom( 'name', 'N' )).xyz[-1]
            Calast = chain.compress(chain.maskFrom( 'name', 'CA' )).xyz[-1]
            Clast = chain.compress(chain.maskFrom( 'name', 'C' )).xyz[-1]
            ## Mejor ir arrastrando los valores de N etc
            
            chain, R = planarize(chain,calcPlane(Nlast,Calast,Clast))

            ## Center chain
            for j in range(len(chain.xyz)):
                chain.xyz[j] = chain.xyz[j] - Nlast
            
            ## Rotar el residuo a su plano original
            
            res , R= planarize(res,[[1,0,0],[0,0,-1]],True,nextres_reg["plane_next"],False)
            
            ## Sumarle la translacion al residuo 
            for j in range(len(res.xyz)):
                res.xyz[j] = res.xyz[j] + nextres_reg["translation"]
            
            ## Pegarlos
            chain = chain.concat(res)
            
            ## Rotar los puntos 
            
            
            ## Y volver a empezar
            r = rnext
            i = i+1
            N = res.xyz[0]
            # sacar tb valor de Ca, C
        return chain
        

    def createStericChain(self):
       ## Lo mismo, pero en lugar de montarla primero crea una lista de modelos
       ## y los va juntando por backtracking
       pass 
    
    
    
    def __str__(self):
        """
        Returns a string with statistics of the database.
        """
        
        mystring = ""
        mystring+= "Residues Database \n =====================\n"
        mystring+= "Total residues picked:"+str(len(self.reslib))+"\n"
        
        mystring+= "\nTotal per residue: \n =====================\n"
        changeline = 0
        for k in MU.allAACodes():
            mystring+=k+":"+str(len(self.index[k]))+"    "
            changeline+=1
            if changeline == 5:
                mystring += "\n"
                changeline = 0
            
        
        return mystring



##############
## Test
##############
import Biskit.test as BT
import Biskit.tools as T
from Biskit.PDBModel import PDBModel
import os

class Test(BT.BiskitTest):
    """ Test cases for residuePicker"""

    def prepare(self):
        pass

    def cleanUp( self ):
        try:
            #~ #os.rmdir(T.testRoot()+"/polysys/residues_db")
            os.system("rm -rf "+T.testRoot()+"/polysys/residues_db")
        except:
            if self.local:
                print "ERROR: Database folder couldn't be removed."
    
    #~ def test_DB_FileCreation(self):
        #~ """Folder & Files creation testing"""
        
        #~ r = residuePicker(T.testRoot()+"/polysys/residues_db/")
        #~ r.extractFromPDB(T.testRoot()+"/polysys/FAKE.pdb")
        #~ r.save()
        
        #~ r2 = residuePicker(T.testRoot()+"/polysys/residues_db/")
        
        #~ ## r and r2 must be the same
        #~ self.assertEqual(r.index,r2.index)
        
    #~ def test_DB_Processing(self):
        #~ """DB Processing test cases"""
        
        #~ r = residuePicker(T.testRoot()+"/polysys/residues_db/")
        #~ ## empty DB
        #~ if self.local:
            #~ print r
  
        #~ r.extractFromPDB(T.testRoot()+"/polysys/FAKE.pdb")
        #~ ## DB with FAKE protein processed
        #~ if self.local:
            #~ print r
       
        #~ len1= len(r.index)
        #~ lendata1= lendepth(r.apvect)
        
        #~ r.extractFromPDB(T.testRoot()+"/polysys/FAKE.pdb")
        #~ ## Just the same DB!!
        #~ if self.local:
            #~ print r
        #~ len2= len(r.index)
        #~ lendata2= lendepth(r.apvect)
        
        #~ self.assertEqual(len1,len2)
        #~ self.assertEqual(lendata1,lendata2)
    
    #~ def test_DB_FileCreation2(self):
        #~ """ Extraction and processing Test"""
        #~ r = residuePicker(T.testRoot()+"/polysys/residues_db/")
        #~ r.extractFromPDB(T.testRoot()+"/polysys/1HUY.pdb")
        #~ if self.local:
            #~ print r
        #~ path = T.testRoot()+"/polysys/residues_db/K/"
        #~ files = os.listdir(path)
        
        #~ f = open(T.testRoot()+"/polysys/tots.pdb","w")
        
        #~ i=1
        #~ for r in files:
            #~ f2 =  open(path+r,"r")
            #~ f.write("MODEL    %4d\n"%(i))
            #~ f.writelines(f2.readlines()[:-1])
            #~ f.write("ENDMDL\n")
            #~ f2.close()
            #~ i=i+1
        #~ f.close
            
    #~ def test_randomres(self):
        #~ """Test random res picker function"""
        #~ r = residuePicker(T.testRoot()+"/polysys/residues_db/")
        #~ r.extractFromPDB(T.testRoot()+"/polysys/FAKE.pdb")
        #~ if self.local:
            #~ print r.randomRes('L')
            #~ print r.randomRes(r.index['L'])
        #~ r.extractFromPDB(T.testRoot()+"/polysys/1HUY.pdb")
        #~ if self.local:
            #~ print r.randomRes('A')
            #~ print r.randomRes(r.index['A'])
            #~ print r.randomRes('A')
            #~ print r.randomRes(r.index['A'])
            #~ print r.randomRes('A')
            #~ print r.randomRes(r.index['A'])
            #~ print r.randomRes('A')
            #~ print r.randomRes(r.index['A'])
            #~ print r.randomRes('A')
            #~ print r.randomRes(r.index['A'])
    
    def test_chaincreation(self):
        """Test chain creation function"""
        r = residuePicker(T.testRoot()+"/polysys/residues_db/")
        p = PDBModel(T.testRoot()+"/polysys/FAKE3.pdb")
        #~ r.verbose = True
        r.extractFromPDB(p)
        chain = r.createChain("DPM",genPoints(p))
        chain.writePdb(T.testRoot()+"/polysys/chain.pdb")
        
if __name__ == '__main__':
    BT.localTest()    




