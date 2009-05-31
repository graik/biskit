from Biskit.PDBModel import PDBModel
from numpy import compress,transpose,where,cross,array,matrix,array
import cPickle 
import os
from math import acos,sqrt,atan2
import emath 
from emath import vectorangle
from quaternion import rotquat, rotmat
import Biskit.molUtils as MU 


class residuePicker:
    
    def __init__(self,databasepath = "./residues_db/"):
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
        
        """
        
        self.dbpath = databasepath
        
        try:
            print "db  ", self.dbpath+'index.db'
            f = open (self.dbpath+'index.db',"r")
            print f.readlines()
            self.index = cPickle.load(f)

        except (cPickle.UnpicklingError, IOError,EOFError):
  
            if( not os.access(self.dbpath, os.F_OK)):
                os.mkdir(self.dbpath)
            f = open (self.dbpath+'index.db',"w")
            self.index = {}
        f.close()
        
        try:
            f = open (self.dbpath+'apvect.db',"r")
            self.apvect = cPickle.load(f)

        except (cPickle.UnpicklingError, IOError,EOFError):
            f = open (self.dbpath+'apvect.db',"w")
            self.apvect = {}
        f.close()
        
        try:
            f = open (self.dbpath+'clusters.db',"r")
            self.clusters = cPickle.load(f)

        except (cPickle.UnpicklingError, IOError,EOFError):
            f = open (self.dbpath+'clusters.db',"w")
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
        try:
            f = open (self.dbpath+'index.db',"w")
            cPickle.dump(self.index,f,cPickle.HIGHEST_PROTOCOL)
            f.close()
        except:
            print "ERROR"
        
        try:
            f = open (self.dbpath+'index.db',"r")
            self.index = cPickle.load(f)
        except (cPickle.UnpicklingError, IOError,EOFError):
            print "ERROR2"
          
        
        for k in self.index:
            print k, len(self.index[k])
        
        f = open (self.dbpath+'apvect.db',"w")
        cPickle.dump(self.apvect,f,cPickle.HIGHEST_PROTOCOL)
        f.close()
        
        f = open (self.dbpath+'clusters.db',"w")
        cPickle.dump(self.clusters,f,cPickle.HIGHEST_PROTOCOL)
        f.close()
        
        print self.clusters

    def extractFromPDB ( self,path = "" ):
        """
        
        """
        ## Clean the structure a bit
        model = PDBModel(path)
        model = model.compress( model.maskProtein())
        
        ## Also remove all hidrogens
        model = model.compress(model.maskHeavy())
        

        model.report()
        for c in range( model.lenChains() ):
            
            chain = model.takeChains( [c] )
            
            residuos = chain.resModels()
            
            for i in range(len(residuos)) :
                print "-"+ residuos[i].sequence() + str(self.last_res[residuos[i].sequence()])
                r = residuos[i].sequence()
                residue = residuos[i]
                myindex = r+str(self.last_res[r])
                
                ## Fetch the next residue
                if i<len(residuos)-1 :
                    nextres = residuos[i+1]
                    rnext = nextres.sequence()
                else:
                    nextres = None
                
                ## Hunt for N & C to reorient the aa
                iNC= where(residue.maskFrom( 'name', ['N','CA','C'] ))
                if nextres != None:
                    iNC2= where(residue.maskFrom( 'name', ['N','CA','C'] ))
                else:
                    iNC2 = None
                    
                ## First purge: it has all the needed atoms?
                if len(residue.xyz) == len(MU.aaAtoms[MU.single2longAA(r)[0]])-1:
                                        
                    (residue,vector) = self.doReorientation(residue,iNC,nextres,iNC2)
 
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
                                    print "Too similar to other ("+ similarto +")"
                                    
                                else:##if not foundsd:
                                    group.append(myindex)
                                    add_it = True
                                    print "Backbone found but different sidechain:adding"
                                    
                    if not foundback:
                        ## We haven't found the backbone so we need to create a new
                        ## group.
                        self.clusters[r].append([myindex])
                        add_it = True
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
                    
                    if nextres  != None :
                        key = (r,rnext)
                        new_entry =  (myindex,vector)
                        if not key in self.apvect.keys():
                            ## if the key isn't we just add the new register
                            self.apvect[key] = [ new_entry ]   
                            
                        else:
                            ## else we have to see if there's something similar inside
                            foundsimilar = False
                            for i in self.apvect[key]:
                                
                                if vectorangle(vector[0],i[1][0]) < 0.045:
                                    foundsimilar = True
                            if foundsimilar == False:
                                self.apvect[key].append(new_entry)
                else:
                    print "Diferent number of atoms."
        self.save()
    
    def doReorientation (self,residue = None, iNC = [0,1,2], nextres = None, iNC2 = [0,1,2]):
        """
        Reorients a residue along the C-Ca axis and C-N axis.
        It's used for having an standard base orientation for all the residues.
        
        @param residue: Residue to be reoriented.
        @type residue: PDBModel
        @param iNC: Atom indexes of the N, CA,C atoms.
        @type iNC: list of int
        
        @return: The residue reoriented.
        @rtype: PDBModel
        """

        if nextres  != None :
            N2 = nextres.xyz[iNC2[0][0]]

        N = residue.xyz[iNC[0][0]]
        if nextres  == None :
            N2 = N
        
        C = residue.xyz[iNC[0][2]]
        vector = N2 - C
        
        ## Center
        residue.xyz = residue.xyz - N
        
        ## Rotation
        C = residue.xyz[iNC[0][2]]
        
        
        
        Ca = residue.xyz[iNC[0][1]]
        nNCNCa = emath.normalized(cross(C,Ca))	
        NC = array([[0.,0.,0.],[0.,0.,0.],nNCNCa])
        NC2 = array([[0.,0.,0.],[0.,0.,0.],[0.,0.,-1.]])	
        q= rotquat(NC,NC2)[2]
        R = transpose(matrix(rotmat(q)))				
        for i in range(len(residue.xyz)):
            residue.xyz[i] = residue.xyz[i] * R
        residue.update()
        
        vector = vector * R
        
        C = residue.xyz[iNC[0][2]]
        Cnorm = emath.normalized(C)
        NC = array([[0.,0.,0.],[0.,0.,0.],Cnorm])
        NC2 = array([[0.,0.,0.],[0.,0.,0.],[0.,1.,0.]])
        q= rotquat(NC,NC2)[2]
        R = transpose(matrix(rotmat(q)))
        for i in range(len(residue.xyz)):
            residue.xyz[i] = residue.xyz[i] * R
        residue.update()
        
        vector = vector * R
        

        return residue,array(vector)
        
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

    
    def returnRandomRes(self, resname = 'A'):
        """
        Returns a random residue loaded from the library of type resname.
        
        @param resname: Single letter code of aa.
        @type resname: string
        
        @return: True if the too residues are very similar in conformation.
        @rtype: PDBModel
        """
        pass
        
        
    def rotateAA (self, paa='',v1 = [0.,0.,0.],v2 = [0.,0.,0.],writeIt = False,name=''):
        
        if not isinstance(paa,PDBModel):
            aa = PDBModel(paa)
        else:
            aa = paa
            
        # we assume it's in the YZ plane ( it must be from DB )
        
        NC = array([[0.,0.,0.],[0.,0.,0.],v1])
        NC2 = array([[0.,0.,0.],[0.,0.,0.],v2])
        q= rotquat(NC,NC2)[2]
        #~ print "q:",q,rotquat(NC,NC2)
        R = transpose(matrix(rotmat(q)))	
        #~ print R,v1,v2
        for i in range(len(aa.xyz)):
            aa.xyz[i] = aa.xyz[i] * R
        aa.update()
        
        if writeIt :
            aa.writePdb(name)
        
        return aa

    def stickAAs( self,a,b , flip ):
        iC= where(a.maskFrom( 'name', ['C'] ))
        #~ print a
        #~ print b
        
        if flip :
            
            b = self.rotateAA(b,[1.,0.,0.],[0.,0.,1.],True,'lol.pdb')
            b = self.rotateAA(b,[0.,0.,1.],[-1.,0.,0.],True,'lol.pdb')
            
            T = array([0,0.56,1.27])
        
        else:
            T = array([0,-0.56,-1.27])
            
        C = a.xyz[iC[0][-1]]
        
        for i in range(len(b.xyz)):
            b.xyz[i] = b.xyz[i] + T + C
        b.update()
        
        
        m = a.concat(b)
        #~ print m
        m.renumberResidues()
        
        m['serial_number'] = range(1,m.lenAtoms()+1)
        m['chain_id'] = m.lenAtoms() * ['A']
        
        return m
        
    def concatenateAAChains (self,a,b) :
        # temptative function for chain concatenation
        # chains may not have TER o OXT
        
        #reorientation of the other molecule in the or. of the other (Nn->Cn to b's N1->C1)
        iCa= where(a.maskFrom( 'name', ['C'] ))
        iNa= where(a.maskFrom( 'name', ['N'] ))
        
        iCb= where(b.maskFrom( 'name', ['C'] ))
        iNb= where(b.maskFrom( 'name', ['N'] ))
        
        #Residues are supposed to have the correct backbone at least
        aCn = a.xyz[iCa[0][-1]] # Cn
        aNn = a.xyz[iNa[0][-1]] # Nn
        
        aCnNn = emath.normalized(aCn - aNn)
        
        bC1 = b.xyz[iCb[0][0]] # C1
        bN1 = b.xyz[iNb[0][0]] # N1
        
        bC1N1 = emath.normalized(bC1 - bN1)
        
        #center
        a.xyz = a.xyz - aNn
        b.xyz = b.xyz - bN1
        #~ print "centered-",bN1,bC1,iCb[0][0],iCa[0][0]
        #Reorient b 
        b = self.rotateAA(b,bC1N1,aCnNn)
        
        c =self.stickAAs(a,b,False)
        #~ print a 
        #~ print b
        #~ print c
        return c


    def calculateBulkRMS(self,aaname = 'A', writeIt = False):
        
           
        if writeIt:
            f = open ("./rmsd_"+aaname,"w")
        
        #load all aa of name aaname
        for i in range(1,total_res[aaname]+1):
            print "./residues_db/"+aaname+"/"+aaname+str(i)
            residues.append(PDBModel("./residues_db/"+aaname+"/"+aaname+str(i)+".pdb"))
            
        
        
        for i in range(total_res[aaname]):
            for j in range(total_res[aaname]):
                if i!= j:
                    if writeIt:
                        f.write(str(i)+" "+str(j)+ " " + str (residues[i].rms(residues[j]))+"\n")
            f.write("\n")
        
        if writeIt:
            f.close()
        
        

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
    """ Test cases for Polyfret"""

    def prepare(self):
        pass

    def cleanUp( self ):
        #~ os.system("rm -rf "+T.testRoot()+"/polysys/residues_db")
        pass
    
    def test_DB_Creation(self):
        """DB Creation test cases"""
        print T.testRoot()+"\\polysys\\residues_db\\"
        r = residuePicker(T.testRoot()+"\\polysys\\residues_db\\")
        
        if self.local:
            print r
            
            
        r.extractFromPDB(T.testRoot()+"/polysys/FAKE.pdb")
        if self.local:
            print r
        #~ print r.apvect
        r.extractFromPDB(T.testRoot()+"/polysys/FAKE.pdb")
        
        if self.local:
            print r
        #~ print r.apvect 
        #~ print vectorangle([-1.11594935, -0.00334113, -0.72240529] ,[-1.11594935, -0.00334113, -0.72240529])

if __name__ == '__main__':
    BT.localTest()    




