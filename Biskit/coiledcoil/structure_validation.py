from validation import lengthFilter,file_stats
from Biskit import PDBModel
import os
from getcoil import createCandidatesFile,dataFileCreation,getBestHit,doHomologyModelling
from math import sqrt
from Biskit.tmalign import TMAlign


def getCoilStructs ( basepath ,candidates_file):
    file = open(candidates_file)
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
        if not contents[0][:-4]+"@"+contents[3] in coils:
            ## pdb@helixId
            coils[contents[0][:-4]+"@"+contents[3]] = [model.takeResidues(range(int(contents[4])-1,int(contents[5])))]
        else:
            coils[contents[0][:-4]+"@"+contents[3]].append( model.takeResidues(range(int(contents[4])-1,int(contents[5]))))
        correspondences[contents[0][:-4]+"@"+contents[4]]= contents[0][:-4]+"@"+contents[3]
    
    structs = {} 
    
    sequences = {}
    chains = "ABCDEFGHI"
    for coil in coils:
        pdb = None
        i = 0
        #~ print "processing coil",coil
        for helix in coils[coil]:
            #~ print i
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

def massiveRMSD(pdb1,pdb2):
    ## Only for CA 
    seq1 = pdb1.sequence()
    seq2 = pdb2.sequence()
    
    pdb1.renumberResidues()
    pdb2.renumberResidues()
    
    #~ print "longitudes",len(seq1),len(seq2)
    #~ print seq1
    #~ print seq2
    if len(seq1) != len(seq2):
        if len(seq1) >= len(seq2):
            a = pdb1
            b = pdb2
        else:
            a = pdb2
            b = pdb1
        
        extra = abs(len(seq1)-len(seq2))
        
        rmsd = []
        #~ print "extra",extra
        for i in range(extra+1):
            a_=a.takeResidues(range(i,i+len(b.sequence())))
            #~ print a_.sequence()
            #~ print b.sequence()
            rmsd.append((a_.rms(b),a.concat(b)))
        #~ print rmsd    
        return min(rmsd)[0],len(b)
    else:
        return pdb1.rms(pdb2),float(len(pdb1.atoms))

def massiveDimerRMSD(pdb1,pdb2):
    ## Only for CA 
    seq1 = pdb1.sequence()
    seq2 = pdb2.sequence()
    
    #~ print "seq1", seq1
    #~ print "seq2", seq2
    
    pdb1.renumberResidues()
    pdb2.renumberResidues()
    
    pdb1_chainA = pdb1.takeChains([0]) 
    pdb1_chainB = pdb1.takeChains([1])
    pdb2_chainA = pdb2.takeChains([0])
    pdb2_chainB = pdb2.takeChains([1])
    
    seq1A = pdb1_chainA.sequence()
    seq1B = pdb1_chainB.sequence()
    seq2A = pdb2_chainA.sequence()
    seq2B = pdb2_chainB.sequence()
    
    changedA = changedB = False
    
    #~ print seq1A, seq1B
    #~ print seq2A, seq2B
    
    if len(seq1A) >= len(seq2A):
        _1a = pdb1_chainA
        _2a = pdb2_chainA
    else:
        changedA = True
        #~ print "changing A"
        _1a = pdb2_chainA
        _2a = pdb1_chainA
    
    if len(seq1B) >= len(seq2B):    
        _1b = pdb1_chainB
        _2b = pdb2_chainB
    else:
        changedB = True
        #~ print "changing B"
        _1b = pdb2_chainB
        _2b = pdb1_chainB

    if len(seq1) != len(seq2):
                
        extraA = abs(len( _1a.sequence()) - len(_2a.sequence()))
        extraB = abs(len( _2a.sequence()) - len(_2b.sequence()))
        
        rmsd = []
        #~ print "First dimer",_1a.sequence(), _2a.sequence()
        #~ print "Second dimer",_1b.sequence(), _2b.sequence()
        
        #~ print "extra",extra
        for i in range(extraA+1):
            s_1a=_1a.takeResidues(range(i,i+len(_2a.sequence())))
            #~ print a_.sequence()
            #~ print b.sequence()
      
            for j in range(extraB+1):
                
                s_1b=_1b.takeResidues(range(j,j+len(_2b.sequence())))
                #~ print s_1a.sequence(), s_1b.sequence()
                #~ print _2a.sequence(), _2b.sequence()
                
                if changedA and changedB:
                    ## nothing happens
                    pass
                else:
                    if changedB:
                        #~ print "reordering Bs"
                        aux = s_1b
                        s_1b = _2b
                        _2b = aux
                    
                    if changedA:
                        #~ print "reordering As"
                        aux = s_1a
                        s_1a = _2a
                        _2a = aux
                
                  
                total1 = s_1a.concat(s_1b)
                total2 = _2a.concat(_2b)
                    
                rmsd.append(total1.rms(total2))
        #~ print rmsd    
        return min(rmsd)
    else:
        return pdb1.rms(pdb2)

def dimerMinRmsd(pdb1,pdb2):
    
    pdb1_chainA = pdb1.takeChains([0]) 
    pdb1_chainB = pdb1.takeChains([1])
    pdb2_chainA = pdb2.takeChains([0])
    pdb2_chainB = pdb2.takeChains([1])
        
    rmsaa,NAA = massiveRMSD(pdb1_chainA,pdb2_chainA)
    rmsab,NAB = massiveRMSD(pdb1_chainA,pdb2_chainB)
    rmsba,NBA = massiveRMSD(pdb1_chainB,pdb2_chainA)
    rmsbb,NBB = massiveRMSD(pdb1_chainB,pdb2_chainB)
    
    dsumAA = NAA*rmsaa
    dsumAB = NAB*rmsab
    dsumBA = NBA*rmsba
    dsumBB = NBB*rmsbb
    
    return sqrt((dsumAA+dsumBB)/(NAA+NBB)),sqrt((dsumAB+dsumBA)/(NAB+NBA)),min(sqrt((dsumAA+dsumBB)/(NAA+NBB)),sqrt((dsumAB+dsumBA)/(NAB+NBA)))
    
    

def validate(data_path,candidates_file):
    pdb_used = []
    
    structs, sequences, correspondences= getCoilStructs(basepath,candidates_file)
   
    data_val_path = data_path+"_new"
    
    rmsd_total = {}
    
    ## Do a massive RMSD
    
    
    for s1 in structs:
        for s2 in structs:
            mdrmsd = massiveDimerRMSD(structs[s1].compress(structs[s1].maskCA()),structs[s2].compress(structs[s2].maskCA()))
            tmaligner = TMAlign( structs[s1], structs[s2] )
            res = tmaligner.run()
            rmsd_total[(s1,s2)]= {'shiftRMSD':mdrmsd,'tmalRMSD':res['rmsd'],'tmalTSCO':res['score']}
    

    rmsd = {}
    for pdb in structs:
        
        ## Extract one from data file and retarget
        print "-",pdb
        deleted = delFromDataFile(data_path,data_val_path,pdb)
        
        if deleted:
            ## Predict structure
            #~ print data_val_path, sequences[pdb]
            study, best = getBestHit ( datafile =data_val_path, method = "Parry",sequences=sequences[pdb])
            #~ print study, best
            if best[1] != None:
                #~ print pdb, best[1][1], sequences[pdb]
                templa_path,result_path = doHomologyModelling ( id = best[1][1] ,candidates_file = candidates_file, sequences = sequences[pdb], study = study)
                ## Get the result
                ## (find the file)
                files = os.listdir(result_path)
                for f in files:
                    if ".pdb" in f and "model_" in f:
                        resultpdb = PDBModel(result_path+"/"+f)
                
                
                ## Make rmsd x.magicFit(p)

                structs[pdb].writePdb(result_path+"/original.pdb")
                myrmsd = structs[pdb].compress(structs[pdb].maskCA()).rms(resultpdb.compress(resultpdb.maskCA()))
                tmaligner = TMAlign( structs[pdb], resultpdb )
                res = tmaligner.run()
                rmsd[pdb] = {'RMSD':myrmsd,'tmalRMSD':res['rmsd'],'tmalTSCO':res['score'],'selected':best[1][1]}
                print "RMSD",rmsd[pdb]
                ## Clean everything
                os.system("rm -rf "+templa_path)
            else:
                print pdb, "was unpredictable with this set"
               
            #~ print best, pdb
            #~ return
        else:
            print pdb, "is not inside the data file"

    
    print "RMSDS TOTAL--------"
    for r in rmsd_total:
        print "%s -  shiftRMSD:%.3f tmalRMSD:%.3f tmalTSCO:%.3f"%(r, rmsd_total[r]['shiftRMSD'],rmsd_total[r]['tmalRMSD'],rmsd_total[r]['tmalTSCO']) 
    
    print "RMSDS--------"
    for r in rmsd:
        print "%s (%s)-  RMSD:%.3f tmalRMSD:%.3f tmalTSCO:%.3f  %.3f"%(r, rmsd[r]['selected'],rmsd[r]['RMSD'],rmsd[r]['tmalRMSD'],rmsd[r]['tmalTSCO'],abs(rmsd[r]['RMSD']-rmsd[r]['tmalRMSD'])) 
    
        

def extractPdbs(data_path):
    
    file = open(data_path)
    lineas = file.readlines()
    file.close()
    lineas = [ l.strip() for l in lineas ]
    pdbs = []
    for l in lineas:
        contents = l.split()
        pdbs.append( contents[1][7:])
    pdbs_s= set(pdbs)
    for p in pdbs_s:
        print "cp ./pdbs/"+p+" ./pdbs_nano/"+p 
basepath = '/home/victor/poly0.5/Biskit/testdata/coiledcoil/pdbs_nano'

#~ extractPdbs(basepath+"/datos_new")
#~ createCandidatesFile(basepath+"/candidates",basepath)
#~ structs, sequences, correspondences= getCoilStructs(basepath)
#~ print structs
#~ print sequences
#~ dataFileCreation( data_file = basepath+"/datos",candidates_file = basepath+"/candidates", target_seq = 'LLLLLLLLLLL', method = "Any",target_type = ("hetero","antiparallel",2))
   
#~ lengthFilter("/home/victor/poly0.5/Biskit/testdata/coiledcoil/candidates","/home/victor/poly0.5/Biskit/testdata/coiledcoil/candidates_filtered",7)
#~ file_stats("/home/victor/poly0.5/Biskit/testdata/coiledcoil/pdbs_small/candidates")
#~ file_stats("/home/victor/poly0.5/Biskit/testdata/coiledcoil/candidates_filtered")
validate(basepath+"/datos",basepath+"/candidates")
#~ delFromDataFile(basepath+"/datos",basepath+"/datos_new","2FHA")