from validation import lengthFilter,file_stats
from Biskit import PDBModel
import os
from getcoil import createCandidatesFile,dataFileCreation,getBestHit,doHomologyModelling

def getCoilStructs ( basepath ):
    file = open(basepath+'/candidates')
    lineas = file.readlines()
    file.close()
    lineas = [ l.strip() for l in lineas ]
    pdbs_path = lineas[0]
    coils = {}
    correspondences = {}    
    for l in lineas[1:]:
        contents = l.split()
        print "Processing..."+contents[0]
        model = PDBModel(pdbs_path+"/"+contents[0])
        ### Do same preprocessing as in socket 
        mask  = model.maskProtein(standard=True)
        model = model.compress( mask )
        model.renumberResidues()
        if not contents[0]+contents[3] in coils:
            ## pdb@helixId
            coils[contents[0][:-4]+"@"+contents[3]] = [model.takeResidues(range(int(contents[4]),int(contents[5])+1))]
        else:
            coils[contents[0][:-4]+"@"+contents[3]].append( model.takeResidues(range(int(contents[4]),int(contents[5])+1)))
        correspondences[contents[0][:-4]+"@"+contents[4]]= contents[0][:-4]+"@"+contents[3]
    
    structs = {} 
    
    sequences = {}
    chains = "ABCDEFGHI"
    for coil in coils:
        pdb = None
        i = 0
        for helix in coils[coil]:
            helix.atoms['chain_id']=[chains[i]]*len(helix)

            i = i+1
            if pdb == None:
                pdb = helix
                sequences[coil]= [helix.sequence()]
            else:
                pdb = pdb.concat(helix)
                sequences[coil].append(helix.sequence())
        
        structs[coil]= pdb
    return structs, sequences,correspondences


def delFromDataFile(data_path,data_val_path,who):
    file = open(data_path)
    lineas = file.readlines()
    file.close()
    lineas = [ l.strip() for l in lineas ]
    
    file2 = open(data_val_path,"w")
    
    deleted = False
    target = []
    for l in lineas:
        contents = l.split()
        #~ print contents[0][0:4],who[0:4]
        
        if contents[0][0:4]!=who[0:4] and contents[0][0:4]!='TARG':
            file2.write(l+"\n")
        
    
        if contents[0][0:4] == who[0:4] and target == []:
            target = contents[2:]
            deleted = True
            
    
    file2.write("TARGET ")
    for c in target:
        file2.write(c+" ")
    
    file2.close()
    
    return deleted
    
    
def validate(data_path,candidates_file):
    pdb_used = []
    
    structs, sequences, correspondences= getCoilStructs(basepath)
   
    data_val_path = data_path+"_new"
    
    for pdb in structs:
        
        ## Extract one from data file and retarget
        print "-",pdb
        deleted = delFromDataFile(data_path,data_val_path,pdb)
        
        if deleted:
            ## Predict structure
            print data_val_path, sequences[pdb]
            study, best = getBestHit ( datafile =data_val_path, method = "Parry",sequences=sequences[pdb])
            print study, best
            if best[0] != None and best[1] != None:
                doHomologyModelling ( id = best[1][1] ,candidates_file = "",  sequences = sequences[pdb],study = study)
                ## Get the result
                
                resultpdb = PDBModel()
                
                ## Make rmsd
                #~ rmsd[pdb] = structs[pdb].rms(result)
            else:
                print pdb, "was unpredictable with this set"
                #~ return 
            print best, pdb
            return
        else:
            print pdb, "is not inside the data file"
        

basepath = '/home/victor/poly0.5/Biskit/testdata/coiledcoil/pdbs'

#~ createCandidatesFile(basepath+"/candidates",basepath)

#~ dataFileCreation( data_file = basepath+"/datos",candidates_file = basepath+"/candidates", target_seq = 'LLLLLLLLLLL', method = "Any",target_type = ("hetero","antiparallel",2))
   
#~ lengthFilter("/home/victor/poly0.5/Biskit/testdata/coiledcoil/candidates","/home/victor/poly0.5/Biskit/testdata/coiledcoil/candidates_filtered",7)
#~ file_stats("/home/victor/poly0.5/Biskit/testdata/coiledcoil/pdbs_small/candidates")
#~ file_stats("/home/victor/poly0.5/Biskit/testdata/coiledcoil/candidates_filtered")
validate(basepath+"/datos",basepath+"/candidates")
#~ delFromDataFile(basepath+"/datos",basepath+"/datos_new","2FHA")