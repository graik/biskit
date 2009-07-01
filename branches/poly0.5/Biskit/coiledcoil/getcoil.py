from coiledcoil import CoiledCoil
from choosecoil import CCStudy
import methods
import os
from socket_coil import SocketCoil
from Biskit import PDBModel
from coiledutils import getHeptad

CHAINS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

def createCandidatesFile (filename = "", dir = ""):
    """
    Creates the candidates file.
    It reads a directory with pdb files containing possible coiled coils.
    For each one, it uses socket to create a 
    Please see 'dataFileCreation' for candidates file format description.
    
    @param path: File name of the candidates file.
    @type path:string
    @param dir: Path of the directory containing pdb template files
    @type dir: string
    
    
    @Candidates file format:
    The first line must be the path to a folder containing the pdb structure files
    of our templates.
    
    Then it's followed by several lines with this syntax:
    
    PDB_FILE_NAME OL SENSE ID chain register_representative
    
    Where:
    
    - PDB_FILE_NAME is a string with the file name of the pdb structure model we refer
    
    - OL is the oligomerization state of the coiled coil
    
    - SENSE can be "parallel" or "antiparallel"
    
    - ID is the id of the coiled coil (just to get the complementary chains)
    
    - chain is the chain to be used as reference for alignment
    
    - register_representative is... the register representative, so a 7 residues sentence
        which is defined by the "abcdefg" register.
    """
    
    file = open(filename,"w")
    
    assert (os.path.exists(dir)), "Error, path ( %s ) doesn't exist"%(path)
        
    pdbs = os.listdir(dir)
    
    assert (pdbs != []),"Error, dir provided ( %s ) is empty"%(path)
    
    failed = []
    
    for pdb in pdbs:
        try:
            print "preparing"+ dir+"/"+pdb
            
            model = PDBModel(dir+"/"+pdb)
            sc = SocketCoil(model)
            sc.debug = True
            sc.run()
            
            assert (sc.result != {})
            
            #~ for cc in sc.result.keys():
                #~ print sc.result[cc]
                
            for cc in sc.result.keys():
                line = ""
                if len(sc.result[cc].ranges) >0:
                    
                    for h in sc.result[cc].ranges:
                        line = pdb +" "+str(sc.result[cc].oligomerization)+" "+sc.result[cc].sense+" "
                        line += sc.result[cc].coilId+" "
                        if len(sc.result[cc].chains[h])>7:
                            line += sc.result[cc].ranges[h][0] + " " +sc.result[cc].ranges[h][1] + " "
                            line += sc.result[cc].chains[h] + " " +getHeptad(sc.result[cc].chains[h],sc.result[cc].registers[h])+" "
                        
                        file.writelines(line+"\n")
            
            del sc 
        except:
            print "failed"
            failed.append(pdb)
    
    file.close()


def dataFileCreation( candidates_file = ""  ,target_seq = "",target_type = ("homodimer","parallel"), method = "All"):
    """
    Parses the candidates file and gets the master data file.
    
    @param candidates_file: Path for the candidates file.
    @type candidates_file: string
    @param target_seq: Target sequence.
    @type target_seq: string
    @param target_type: Type of the coiled coil (if dimer):'parallel' or 'antiparallel'
    @type target_type: 'parallel' or 'antiparallel'
    @param method: 'All','Any' or just the name of the method as defined in 
        methods.'METHODS' .
    @type method: string
    """
    
    try:
        file = open(T.dataRoot() + '/coiledcoil/'+table,"r")
        lineas = file.readlines()
    except IOError, msg:
        raise BiskitError('Cannot open score file %s.\n Error: %s' \
                          % (db, msg ) )
    file.close()
    
    assert (not(target_type[0]=="homodimer" or target_type[0]=="heterodimer")),"Options for target type are 'heterodimer' or 'homodimer'." 
    assert (target_seq == ""),"You have to define a target sequence."
    
    lineas = [ l.strip() for l in lineas ]
    basepath = lineas[0]
    
    
    if not try_all:
        if target_type[1] == "parallel":
            ccdb = db or T.dataRoot() + '/coiledcoil/SOCKET_par_norm'
        elif:
            ccdb = db or T.dataRoot() + '/coiledcoil/SOCKET_antipar_norm'
        else:
            raise BiskitError("Options for target type are 'parallel' or 'antiparallel'." )
        cc = CoiledCoil(ccdb)
        
    file = open(basepath+"/_data","w")
    for l in lineas[1:]:
        line = l.split()
        new_line = line[0][:-3]+" struct:"+line[0]
        if try_all:
            for k in METHODS.keys():
                new_line = new_line +" "+getRegisterByMethod(seq,k)
        else:
            new_line = new_line + " Default:" + cc.findHeptads(seq)['best']
        file.writeline(new_line)
    file.close()


#~ def getCoilStructure ( type = ("parallel","homodimer"), db =""):
    #~ """
    #~ Function for coiled coil prediction.
    #~ """
    #~ assert (candidates_file != ""), " A candidates data file is needed for this function to work"
    
    #~ if type[0] == "parallel":
        #~ ccdb = db or T.dataRoot() + '/coiledcoil/SOCKET_par_norm'
    #~ elif:
        #~ ccdb = db or T.dataRoot() + '/coiledcoil/SOCKET_antipar_norm'
    #~ else:
        #~ ## ERROR
        #~ pass
        
    
    #~ study = CCStudy(cc)
    
    #~ ## Data file creation (candidates parsing)

##############
## Test
##############
import Biskit.test as BT
import Biskit.tools as T
from Biskit import PDBModel

class Test(BT.BiskitTest):
    """ Test cases for Coiled Coil Utils"""

    def prepare(self):
        pass

    def cleanUp( self ):
        pass
    def test_candidates(self):
        """ testing of candidates file creation"""
        createCandidatesFile(T.testRoot()+"/coiledcoil/candidates",T.testRoot()+"/coiledcoil/pdbs_small")
        
        
if __name__ == '__main__':
    BT.localTest()    


