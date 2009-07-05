from coiledcoil import CoiledCoil
from choosecoil import CCStudy
import methods
import os
from socket_coil import SocketCoil
from Biskit import PDBModel
from coiledutils import getHeptad
from Biskit.Mod import Modeller
from alignment import PirAlignment
CHAINS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

def createCandidatesFile (filename = "", dir = "", symetric_offset = 0):
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
    
    assert (pdbs != []),"Error, dir. provided ( %s ) is empty"%(path)
    
    
    file.writelines(dir+"\n")
    
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


def dataFileCreation( data_file = "",candidates_file = ""  ,target_seq = "",target_type = ("homo","parallel",2), method = "All"):
    """
    Parses the candidates file and gets the master data file.
    See choosecoil::parseData for more info about the format.
    
    @param dir_save: Absolute path of the file to be saved. If "" then it's saved at the same
        folder used by createCandidatesFile to create the candidates file.
    @type dir_save: string
    @param candidates_file: Path for the candidates file.
    @type candidates_file: string
    @param target_seq: Target sequence.
    @type target_seq: string
    @param target_type: Type of the coiled coil (if dimer):'parallel' or 'antiparallel', then
        'homo'/'hetero' depending in their corresponding partners. The third component
        is the oligomerization state.
    @type target_type: tuple{string,string,int}
    @param method: 'All','Any' or just the name of the method as defined in 
        methods.'METHODS' .
    @type method: string
    """
    
    assert (target_type[0]=="homo" or target_type[0]=="hetero"),"Options for target type are 'heterodimer' or 'homodimer'." 
    assert (target_type[1]=="parallel" or target_type[1]=="antiparallel"),"Options for target type are 'parallel' or 'antiparallel'."
    assert (target_seq != ""),"You have to define a target sequence."
    assert (method in methods.sources + ["All","Any"]),"No valid method ( %s )"%(method)
    
    file = open(candidates_file,"r")
    lineas = file.readlines()
    file.close()
    
    lineas = [ l.strip() for l in lineas ]
    basepath = lineas[0]
    
    if data_file == "":
        file = open(basepath+"/_data","w")
    else:
        file = open(data_file,"w")
        
    for l in lineas[1:]:
        
        contents = l.split()
        
        ## Filtering
        if len(contents[6])>=len(target_seq) and contents[2] == target_type[1]  and int(contents[1]) == target_type[2]:
              
            new_line = contents[0][:-4]+"@"+contents[4]+" struct:"+contents[0]
            new_line += " seq:"+contents[6]+" Socket:"+contents[7]
            
            if method == 'All':
                for k in methods.sources:
                    new_line += " "+methods.getRegisterByMethod(contents[6],k)
            elif method != 'Any':
                new_line +=  " "+methods.getRegisterByMethod(contents[6],method) 
            
            file.writelines(new_line+"\n")
    
    new_line = "TARGET seq:"+target_seq+" "
    
    if method == 'All':
        for k in methods.sources:
            new_line += " "+methods.getRegisterByMethod(target_seq,k)
    elif method != 'Any':
        new_line +=  " "+methods.getRegisterByMethod(target_seq,method) 
    
    file.writelines(new_line+"\n")

    file.close()


def getBestHit ( type = ("parallel","homodimer"), datafile ="",method = "",sequences=[]):
    """
    Function for coiled coil prediction.
    """
    assert (datafile != ""), " A candidates data file is needed for this function to work"
    
    
    cc = methods.getCoilByMethod(method)
    study = CCStudy(data = datafile, cc = cc)
    study.doStudy()
    mystudy= study.chooseBest()
    
    return mystudy,best
    
def doHomologyModelling ( id = "" ,pos = 0, candidates_file = "",  sequences = [],study = None):
    """
    Retrieves the structure predicted by homology modelling.
    
    @param id: 
    @rtype id:
    @param pos: 
    @rtype pos:
    @param candidates_file: Path for the candidates file.
    @type candidates_file: string
    @param sequences: List of the sequences forming part of the coiled coil.
        First sequence in the list must be 'TARGET' sequence.
    @type sequences: List{string}
    """
    
    ## Retrieve complete data from candidates file
    file = open(candidates_file,"r")
    lineas = file.readlines()
    file.close()
    
    lineas = [ l.strip() for l in lineas ]
    basepath = lineas[0]
    
    ## Filter
    
    name = id.split('@')[0]
    first = id.split('@')[1]
    candidates = []
    coilid = ""

    for l in lineas[1:]:
        contents = l.split()
        if contents[0] == name+".pdb" and contents[4] == first:
            coilid = contents[3]
    print coilid
    for l in lineas[1:]:
        contents = l.split()
        if contents[0] == name+".pdb" and contents[3] == coilid:
            candidates.append(l)
    print candidates    
    model = PDBModel(lineas[0]+"/"+name+".pdb")
    
    ### Do same filtering as in socket 
    mask  = model.maskProtein(standard=True)
    model = model.compress(mask )
    model.renumberResidues()
    
    coils = []
    
    templatepath = lineas[0]+"/templates"
    try:
        os.mkdir(templatepath)
    except:
        pass ## the folder is already created
        
    i = 0;
    
    for c in candidates:
        contents = c.split()
        coil = model.takeResidues(range(int(contents[4]),int(contents[5])+1))
        coil.renumberResidues()
        coil.atoms['serial_number'] = range(1,len(coil.atoms['serial_number'])+1)
        coil.atoms['chain_id'] = [CHAINS[i]]*len(coil.atoms['serial_number'])
        coils.append(coil)
        coil.writePdb(templatepath+"/"+CHAINS[i]+".pdb")
        i = i+1
    
    chains = []
    positions = []
    
    ## alignment.pir creation
    al = PirAlignment(chains,positions)
    
    
    chains =[sequences[0]]
    positions = [best[1][0][1]]
    
    ## Use the same aligner
    aligner = mystudy.reg_alignments[id]
    cc = mystudy.mycc
    
    i=1
    data={}
    for s in sequences[1:]:
        data[CHAINS[i]] = {}
        data[CHAINS[i]]["seq"] = coils[i].sequence()
        data[CHAINS[i]]["Socket"] = cc.findHeptads(data[CHAINS[i]]["seq"])[0]
        i=i+1
        
    ## target.fasta creation
    al.writeFasta(templatepath+"/target.fasta")
    
    mod = Modeller(projectFolder = templatepath, resultFolder=templatepath, 
                  fasta_target=templatepath+"/target.fasta", template_folder=templatepath, f_pir=templatepath+"/alignment.pir",
                  starting_model=1, ending_model=i)


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
    
    #~ def test_candidates(self):
        #~ """ testing of candidates file creation"""
        #~ createCandidatesFile(T.testRoot()+"/coiledcoil/candidates_micro",T.testRoot()+"/coiledcoil/pdbs_micro")
    
    def test_datafilecretion(self):
        """ testing of data file creation"""
        dataFileCreation(T.testRoot()+"/coiledcoil/data",T.testRoot()+"/coiledcoil/candidates_micro",target_seq = "ASDFGHVKKLRER")    
        best = getBestHit( datafile = T.testRoot()+"/coiledcoil/data")
        doHomologyModelling ( candidates_file = T.testRoot()+"/coiledcoil/candidates_micro")

if __name__ == '__main__':
    BT.localTest()    


