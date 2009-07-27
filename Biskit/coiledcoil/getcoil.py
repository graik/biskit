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
    
    PDB_FILE_NAME OL SENSE ID chain register_representative HOMO BIG
    
    Where:
    
    - PDB_FILE_NAME is a string with the file name of the pdb structure model we refer
    
    - OL is the oligomerization state of the coiled coil
    
    - SENSE can be "parallel" or "antiparallel"
    
    - ID is the id of the coiled coil (just to get the complementary chains)
    
    - chain is the chain to be used as reference for alignment
    
    - register_representative is... the register representative, so a 7 residues sentence
        which is defined by the "abcdefg" register.
    
    - HOMO can be "homo" or "hetero", defining if chains are all equal or not
    
    """
    print "Starting candidates file generation"
    
    file = open(filename,"w")
    
    assert (os.path.exists(dir)), "Error, path ( %s ) doesn't exist"%(path)
        
    pdbs = os.listdir(dir)
    
    assert (pdbs != []),"Error, dir. provided ( %s ) is empty"%(path)
    
    
    file.writelines(dir+"\n")
    
    failed = []
    total = len(pdbs)
    done = 0
    no_coils = 0
    for pdb in pdbs:
        try:
            print "preparing"+ dir+"/"+pdb
            
            model = PDBModel(dir+"/"+pdb)
            sc = SocketCoil(model)
            sc.debug = True
            sc.run()
            
            if (sc.result == {}):
                no_coils += 1
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
                            line += sc.result[cc].homo 
                        
                            file.writelines(line+"\n")
            
            del sc 
        except:
            print "failed"
            failed.append(pdb)
        done = done+1
        if done%10 ==0 :
            print "Parsed:",done,"Total:",total,"Failed:",len(failed),"( no coils:",no_coils,")"
    
    print "Parsed:",done,"Total:",total,"Failed:",len(failed),"( no coils:",no_coils,")"
       
    file.close()


def dataFileCreation( data_file = "",candidates_file = ""  ,target_seq = "",target_type = ("hetero","parallel",2), method = "All", use_big = 0):
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
    @param use_big: Uses a longer chain than the one described to do heptad searching. A longer
        sequence can help methods to get better hits. If it's set to something higher than 0 then 
        the sequence minimun length becomes 'use_big'.
        A sequence of length lesser than 28 can't be used by paircoil, so forcing a sequence to be like
        this helps to gather more data.
    @type use_big: int
    
    
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
        
        ## Filtering:
        ## Catch all sequences of length bigger or equal than the target one and forming part of a
        ## coiled coil with the same sense (parallel/antiparallel), homology of the chains and 
        ## oligomerization state.
        
        if len(contents[6])>=len(target_seq) and contents[2] == target_type[1]  and int(contents[1]) == target_type[2]\
            and contents[8] == target_type[0]:
              
            new_line = contents[0][:-4]+"@"+contents[4]+" struct:"+contents[0]
            new_line += " seq:"+contents[6]+" Socket:"+contents[7]
            print contents[0]
            if use_big>0:
                
                sequence = _expandSequence(basepath+"/"+contents[0],int(contents[4]),int(contents[5]),use_big)
                realsec = contents[6]
                if method == 'All':
                    for k in methods.sources:
                        reg = methods.getRegisterByMethod(sequence,k)
                        basereg = reg
                        if reg != "":
                            reg = reg.split(":")[1]
                            #~ print "Initial for all",reg,k
                            ## get first position
                            pos = sequence.find(reg)
                            while pos > 0:
                                pos = pos -7
                            
                            #~ print sequence[pos:pos+7]
                            pos = pos +7
                            while not reg in realsec:
                                reg = sequence[pos:pos+7]
                                #~ print reg
                                pos=pos+7
                            
                            reg = k+":"+reg
                            #~ print "Final",reg
                        new_line += " "+reg
                elif method != 'Any':
                    reg = methods.getRegisterByMethod(sequence,k)
                    if reg != "":
                        reg = reg.split(":")[1]
                        
                        ## get first position
                        pos = sequence.find(reg)
                        while pos > 0:
                            pos = pos -7
                        pos = pos +7
                        while not reg in realsec:
                            reg = sequence[pos:pos+7]
                            pos=pos+7
                        
                        reg = k+":"+reg
                    new_line += " "+reg
            
            else: 
                
                sequence = contents[6]
                
                if method == 'All':
                    for k in methods.sources:
                        new_line += " "+methods.getRegisterByMethod(sequence,k)
                elif method != 'Any':
                    new_line +=  " "+methods.getRegisterByMethod(sequence,method) 
            
            file.writelines(new_line+"\n")
    
    new_line = "TARGET seq:"+target_seq+" "
    
    if method == 'All':
        for k in methods.sources:
            new_line += " "+methods.getRegisterByMethod(target_seq,k)
    elif method != 'Any':
        new_line +=  " "+methods.getRegisterByMethod(target_seq,method) 
    
    file.writelines(new_line+"\n")

    file.close()

def _expandSequence(pdb,start,end, extension):
    #~ print "PDB",pdb
    model = PDBModel(pdb)
    ### Do same preprocessing as in socket 
    mask  = model.maskProtein(standard=True)
    model = model.compress(mask )
    model.renumberResidues()
    
    add = 0
    if (end-start +1)<extension:
        add = int((extension - (end-start +1)) / 2 )
    
    if start - add <0 : 
        start = 0
    else:
        start = start-add
    
    if end + add >len(model.sequence()): 
        end = len(model.sequence())
    else:
        end = end+add
    
    coil = model.takeResidues(range(start,end))
    
    #~ print coil.sequence()
    return coil.sequence()

def getBestHit ( datafile ="",method = "",sequences=[]):
    """
    Function for coiled coil prediction.
    """
    assert (datafile != ""), " A candidates data file is needed for this function to work"
    
    
    cc = methods.getCoilByMethod(method)
    study = CCStudy(data = datafile, cc = cc)
    s,a = study.doStudy(True)
    print s
    best= study.chooseBest()
    #~ print best
    return study,best
    
def doHomologyModelling ( id = "" ,candidates_file = "",  sequences = [],study = None):
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
    #~ print coilid
    for l in lineas[1:]:
        contents = l.split()
        if contents[0] == name+".pdb" and contents[3] == coilid:
            candidates.append(l)

    model = PDBModel(lineas[0]+"/"+name+".pdb")
    
    ### Do same preprocessing as in socket 
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
    heptads =[]
    for c in candidates:
        contents = c.split()
        ## WARNING
        ## Is socket using a residue range starting from 0 even if the pdb starts
        ## somewhere else? As the preprocessing makes residue range start on 1, range
        ## goes from x-1 to y instead of xto y+1 as would be logical. Errors can appear
        ## here later.
        coil = model.takeResidues(range(int(contents[4])-1,int(contents[5])))
        coil.renumberResidues()
        coil.atoms['serial_number'] = range(1,len(coil.atoms['serial_number'])+1)
        coil.atoms['chain_id'] = [CHAINS[i]]*len(coil.atoms['serial_number'])
        coils.append(coil)
        heptads.append(contents[7])
        coil.writePdb(templatepath+"/"+CHAINS[i]+".pdb")
        #~ print "*", coil.sequence()

        i = i+1
    
    
    ## There must be the same number of coiled coils and sequences
    assert(len(coils) == len(sequences)),"You must provide as sequences as there are in the template coiled coil."
    
    ## Get pairs coil sequence / target sequence
    chains =[]
    positions = []
    
    ## Use the same aligner
    aligner = study.alignments[id]
    cc = study.mycc
    
    i=0
    data={}
    for c in coils:
        data[CHAINS[i]] = {}
        data[CHAINS[i]]["seq"] = coils[i].sequence()
        data[CHAINS[i]]["Socket"] = heptads[i]
        i = i+1
   
    
    #### ATENCION A LOS HEPTADS
    #~ print data        
    i=0
    template = None
    while len(sequences) > 0:
        results = []
        for s in sequences:
            data["TARGET"] = {"seq":s,"Default":cc.findHeptads(s)['best']}
            
            #~ print cc.findHeptads(s)['best'],data["TARGET"]
            study.data =  data
            study.doStudy(True)
            res = study.chooseBest()
            results.append((res[1],s))
            del data["TARGET"]
            i=i+1
        best = max(results)
        
        #~ print best
        chains.append((coils[CHAINS.find(best[0][1])].sequence(),best[1]))
        
        if template == None:
            template = coils[CHAINS.find(best[0][1])].clone()
        else:
            template= template.concat(coils[CHAINS.find(best[0][1])])
            
        
        positions.append(best[0][0][1])
        print positions
        
        sequences.remove(best[1])
        del data[best[0][1]]
    
    
    template.renumberResidues()
    template.atoms['serial_number'] = range(1,len(template.atoms['serial_number'])+1)
    template.atoms['chain_id'] = [' ']*len(template.atoms['serial_number'])
    template.writePdb(templatepath+"/template.pdb")
        
        
    print chains    
    ## alignment.pir creation
    al = PirAlignment(chains,positions)
    
    ## alignment.pir creation
    al.writePir(templatepath+"/alignment.pir","template")
    
    ## target.fasta creation
    al.writeFasta(templatepath+"/target.fasta","target")
    
    mod = Modeller( projectFolder = templatepath,
                    resultFolder=templatepath+"/results", 
                    fasta_target=templatepath+"/target.fasta", 
                    template_folder=templatepath, 
                    f_pir=templatepath+"/alignment.pir",
                    starting_model=1, ending_model=1)

    mod.run()
    
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
        #~ createCandidatesFile(T.testRoot()+"/coiledcoil/candidates",T.testRoot()+"/coiledcoil/pdbs")
        #~ createCandidatesFile("coil","/home/victor/synplexity/synplexity/structures/hypercycle/templates/coiledcoildimantipar")
        createCandidatesFile("coil","/home/victor/synplexity/synplexity/structures/hypercycle/templates/coils")
    
    #~ def test_datafilecretion(self):
        #~ """ testing of data file creation"""
        #~ dataFileCreation(T.testRoot()+"/coiledcoil/data",T.testRoot()+"/coiledcoil/candidates",target_type = ("homo","parallel",2),target_seq = "LEIRAAFLRRRNTALRTRVAELRQRVQRLRNIVSQYETRYGPL",use_big=32)    
        #~ study,best = getBestHit( datafile = T.testRoot()+"/coiledcoil/data")
        #~ print best
        #~ doHomologyModelling ( "2Z5H_clean@52", T.testRoot()+"/coiledcoil/candidates",["LEIRAAFLRRRNTALRTRVAELRQRVQRLRNIVSQYETRYGPL","LEIRAAFLRRRNTALRTRVAELRQRVQRLRNIVSQYETRYGPL"],study)

if __name__ == '__main__':
    BT.localTest()    


