import Biskit.tools as T
from coiledutils import getRegister, getAllHeptads
import random
from getcoil import dataFileCreation
from choosecoil import CCStudy
from Biskit.molUtils import allAA,single2longAA


""" 
Collection of functions for validation and table generation.
"""


def file_stats(path):
    
    """
    It extracts some statistics from one candidates file (see 
    getcoil::createCandidatesFile  for more info about the format). Then
    it prints the statistics.
    
    @param path: Absolute or relative path for the candidates file    
    @type path: string
    """
    
    file = open(path)
    lineas = file.readlines()
    file.close()
    lineas = [ l.strip() for l in lineas ]
    
    pdbs = {}
    coil_oligo = {'2':0,'3':0,'4':0,'5':0,'6':0,'7':0,'8':0}
    parallel_homo = 0
    parallel_hetero = 0
    antip_homo = 0 
    antip_hetero = 0
    seq_length = 0
    
    mypdb = {}

    for l in lineas[1:]:
        
        contents = l.split()

        if len(contents) > 4:
            if contents[0] in pdbs:
                pdbs[contents[0]] += 1./int(contents[1])
            else:    
                pdbs[contents[0]] = 1./int(contents[1])
            
            coil_oligo[contents[1]] += 1./int(contents[1])
            
            seq_length += (1./int(contents[1])) * len(contents[6])
            
            if contents[2] == 'parallel' and contents[8] == 'homo' :
                parallel_homo += 1./int(contents[1])
            elif contents[2] == 'parallel' and contents[8] == 'hetero':
                parallel_hetero += 1./int(contents[1])
            elif contents[2] == 'antiparallel' and contents[8] == 'homo':
                antip_homo += 1./int(contents[1])
            else:
                antip_hetero += 1./int(contents[1])
            
    
    ## Coiled Coils per file
    total_files = len(pdbs.keys())
    total_coils = 0
    for pdb in pdbs:
        total_coils += pdbs[pdb]
    total_coils2 = 0   
    for c in coil_oligo:
        total_coils2 += coil_oligo[c]
    
    print "--------------------"
    print "Stats for ",path
    print "Files parsed from:",lineas[0]
    print "Total Coils (1):",total_coils
    print "Total Coils (2):",total_coils2
    print "Total Files:",total_files
    print "Coils per file:",total_coils / float(total_files)
    print "Oligomerization:",coil_oligo
    print "Mean sequence length:", seq_length / total_coils
    print "P H:",parallel_homo
    print "P Het:",parallel_hetero
    print "A H:",antip_homo
    print "A Het:",antip_hetero 
    
    #~ for p in pdbs:
        #~ if pdbs[p] - 1 -int(pdbs[p]) >0.1:
            #~ print p, pdbs[p]-1 , int(pdbs[p]), pdbs[p] - 1-int(pdbs[p])
        #~ print p, pdbs[p]
        

def splitInFiles(path,oligo = '2'):
    """
    Generates 4 files corresponding to the different combinations between parallel/antiparallel and 
    homo/hetero, filtering for the oligomerization number defined in parameter 'oligo'.
    """
    file = open(path)
    lineas = file.readlines()
    file.close()
    lineas = [ l.strip() for l in lineas ]
    
    file_par_hom = open("homodimeric_parallel_can","w")
    file_anti_hom = open("homodimeric_antiparallel_can","w")
    file_par_het = open("heterodimeric_parallel_can","w")
    file_anti_het = open("heterodimeric_antiparallel_can","w")
    
    file_par_hom.write(lineas[0]+"\n")
    file_anti_hom.write(lineas[0]+"\n")
    file_par_het.write(lineas[0]+"\n")
    file_anti_het.write(lineas[0]+"\n")
    
    for l in lineas[1:]:
        contents = l.split()
        if len(contents) > 4 and contents[1] ==oligo:
            
            if contents[2] == 'parallel' and contents[8] == 'homo' :
                file_par_hom.write(l+"\n")
            elif contents[2] == 'parallel' and contents[8] == 'hetero':
                file_par_het.write(l+"\n")
            elif contents[2] == 'antiparallel' and contents[8] == 'homo':
                file_anti_hom.write(l+"\n")
            else:
               file_anti_het.write(l+"\n")
    
    
    file_par_hom.close()            
    file_anti_hom.close()            
    file_par_het.close()            
    file_anti_het.close()            





def gen_table(path , exclusion_list=[]):
    """
    Generates a probability table using candidates in one file.
    
    A probability table is a table with one line per aminoacid, so:
    
    - Each row starts with the aa. code and then it has 7 float values representing
    the probability of this aa.to be in such position (from a(1)-g(7))
    
    - Each colum has the probabilities of all the aa.s to be in this position.
    
    So for theoretical aa. Documentin (DMT) and Pythonin (PTH) table would be:
    
    (line containing 'abcdefg' is not present in the table)
    
        (a   b   c   d   e   f)
    DMT 0.9 0.0 0.0 0.9 0.5 0.2
    PTH 0.1 0.2 0.3 0.3 0.2 0.1
    
    
    @param path: Path for the candidates file.
    @type path: string
    @param exclusion_list: List of sequence ID's which are not going to be used.
    @type exclusion_list: List{string}
    """
    file = open(path)
    lineas = file.readlines()
    file.close()
    lineas = [ l.strip() for l in lineas ]
    
    table ={}
    for aa in allAA():
        table[aa] = {'a':0,'b':0,'c':0,'d':0,'e':0,'f':0,'g':0}
    
    
    all_heptads = []
    for l in lineas[1:]:
        contents = l.split()
        
        if len(contents[6]) >= 13 and len(contents[7])==7 and not contents[0] in exclusion_list:
            register = getRegister(contents[7],contents[6])
            heptads = getAllHeptads(contents[7],register)
            for h in heptads:
                all_heptads.append(h)
                
    heptads_no_rep = set(all_heptads)
    
    reg = "abcdefg"
    for h in heptads_no_rep:
        for i in range(7):
            table[h[i]][reg[i]] +=1
    
    return table




def normalizeTable(table):
    """
    Normalizes each colum of a table to sum 1. 
    So for the table from 'gen_table' doc:
    
        (a   b   c    d   e   f)
    DMT 0.9 0.0 0.0 0.75 0.7 0.66
    PTH 0.1 1.0 1.0 0.25 0.3 0.33
         1.  1.  1.  1.   1.  1.
    
    """
    reg = 'abcdefg'
    for r in reg:
        total = 0
        for aa in table:
            total += table[aa][r]
        
        for aa in table:
            table[aa][r] = table[aa][r]/float(total)  or 0.

    return table        
    
    
def writeTable(table,path):
    """
    Writes a probability table 'table' to the file "path".
    """
    file = open(path,"w")
    for aa in table:
        line = single2longAA(aa)[0]
        
        for r in "abcdefg":
            line += " %.4f"%table[aa][r]
        
        file.write(line+"\n")
    file.close()


def generateValidationGroups(path,train_percent):
    """
    Generates a random list of Ids with length equals to the 'train_percent' percent of the number
    of sequences in the candidates file.
    
    
    """
    
    file = open(path)
    lineas = file.readlines()
    file.close()
    lineas = [ l.strip() for l in lineas ]
    
    pdbs = []
    for l in lineas[1:]:
        contents = l.split()
        if not contents[0] in pdbs:
            pdbs.append(contents[0])
    print len(pdbs)
    
    random.shuffle(pdbs)

    return pdbs[0:int(train_percent*len(pdbs))],pdbs[int(train_percent*len(pdbs)):]


def lengthFilter(path1,path2,length):
    """
    Writes in 'path2' the filtered candidates file in 'path1' which consists only in 
    sequences greater than 'lenght'(int).
    """
    file2 = open(path2,"w")
    
    file1 = open(path1)
    lineas = file1.readlines()
    file1.close()
    lineas = [ l.strip() for l in lineas ]
    file2.write(lineas[0]+"\n")
    pdbs = []
    for l in lineas[1:]:
        contents = l.split()
        #~ len(contents) > 4 and 
        if len(contents[6])>=length:
             file2.write(l+"\n")
    file2.close()


def incompleteHeptadFilter( path1,path2):
    """
    Writes in 'path2' the filtered candidates file in 'path1' which consists only in 
    sequences where the representative heptad length equals 7.
    """
    file2 = open(path2,"w")
    
    file1 = open(path1)
    lineas = file1.readlines()
    file1.close()
    lineas = [ l.strip() for l in lineas ]
    file2.write(lineas[0]+"\n")
    pdbs = []
    for l in lineas[1:]:
        contents = l.split()
        if len(contents) > 4 and len(contents[7])==7:
             file2.write(l+"\n")
    file2.close()
    


    
def writeValidationCandidates(path1,path2,val=[]):
    """
    Writes a candidates file filtering the sequences in 'path1' by sequences ID in 'val'
    and then writes them to 'path2'.
    """
    file2 = open(path2,"w")
    
    file1 = open(path1)
    lineas = file1.readlines()
    file1.close()
    lineas = [ l.strip() for l in lineas ]
    file2.write(lineas[0]+"\n")
    pdbs = []
    for l in lineas[1:]:
        contents = l.split()
        if contents[0] in val:
            file2.write(l+"\n")
    file2.close()
    
def openTable(path):
    """
    Loads the table in 'path' and returns it.
    """
    file1 = open(path)
    lineas = file1.readlines()
    file1.close()
    lineas = [ l.strip() for l in lineas ]
    r = {}
    for l in lineas:
        aux = l.split()
        values = [ float( n ) for n in aux[1:] ]
        if len(aux[0]) != 3 :
            r[aux[0]] = {'a':values[0],'b':values[1],'c':values[2],'d':values[3],'e':values[4],'f':values[5],'g':values[6]}
        
        else:
            r[single2longAA(aux[0])[0]] = {'a':values[0],'b':values[1],'c':values[2],'d':values[3],'e':values[4],'f':values[5],'g':values[6]}
        
    return r    

def errorFilter(path1,path2):
    """
    Writes in 'path2' the filtered candidates file in 'path1' which consists only in 
    sequences where all the coils are completed (it means, dimers have 2 sequences, trimers 3 and 
    so on... )
    """
    file2 = open(path2,"w")
    
    file1 = open(path1)
    lineas = file1.readlines()
    file1.close()
    lineas = [ l.strip() for l in lineas ]
    file2.write(lineas[0]+"\n")
    pdbs = []
    for l in lineas[1:]:
        
        contents = l.split()
        
        
        if len(contents) > 4 and len(contents[7])==7:
             file2.write(l+"\n")
    file2.close()
    
    
    
#~ import os
## All files
file_stats('coils.db')   

## Avoid incomplete heptads
## 
#~ incompleteHeptadFilter("coils.db","coils_nh.db")
#~ lengthFilter('coils.db','coils_nh.db',28)
#~ file_stats('coils_nh.db')  
#~ splitInFiles('coils.db',oligo = '3') 
#~ file_stats('homodimeric_parallel_can')
#~ lengthFilter('homodimeric_parallel_can','homodimeric_parallel_can20',40)
#~ file_stats('homodimeric_parallel_can20')
#~ file_stats('homodimeric_antiparallel_can')
#~ lengthFilter('homodimeric_antiparallel_can','homodimeric_antiparallel_can20',40)
#~ file_stats('homodimeric_antiparallel_can20')
#~ os.system("cp coils_nh.db all")
#~ splitInFiles('all') 
#~ splitInFiles('coils.db','coils_filteredblah.db',14)
## length filter for candidates
#~ lengthFilter('homodimeric_parallel_val','homodimeric_parallel_val14',14)
#~ lengthFilter('homodimeric_parallel_val','homodimeric_parallel_val14',14)
#~ lengthFilter('homodimeric_parallel_val','homodimeric_parallel_val14',14)
#~ lengthFilter('homodimeric_parallel_val','homodimeric_parallel_val14',14)
#~ lengthFilter('all','all_val4',14)


#~ ## normalize not normalized tables
#~ table= openTable(T.dataRoot()+'/coiledcoil/SOCKET_antipar_norm')
#~ table= openTable(T.dataRoot()+'/coiledcoil/SOCKET_par_norm')
#~ writeTable(normalizeTable(openTable(T.dataRoot()+'/coiledcoil/SOCKET_par_norm')),T.dataRoot()+'/coiledcoil/SOCKET_par_norm_prob')
#~ writeTable(normalizeTable(openTable(T.dataRoot()+'/coiledcoil/SOCKET_antipar_norm')),T.dataRoot()+'/coiledcoil/SOCKET_antipar_norm_prob')



#~ t,v=generateValidationGroups('homodimeric_parallel_can',0.6)
#~ print v
#~ writeTable(normalizeTable(gen_table('homodimeric_parallel_can',v)),T.dataRoot()+'/coiledcoil/homodimeric_parallel')
#~ t,v=generateValidationGroups('homodimeric_antiparallel_can',0.6)
#~ writeTable(normalizeTable(gen_table('homodimeric_antiparallel_can',v)),T.dataRoot()+'/coiledcoil/homodimeric_antiparallel')
#~ t,v=generateValidationGroups('heterodimeric_parallel_can',0.6)
#~ writeTable(normalizeTable(gen_table('heterodimeric_parallel_can',v)),T.dataRoot()+'/coiledcoil/heterodimeric_parallel')
#~ t,v=generateValidationGroups('heterodimeric_antiparallel_can',0.6)
#~ writeTable(normalizeTable(gen_table('heterodimeric_antiparallel_can',v)),T.dataRoot()+'/coiledcoil/heterodimeric_antiparallel')
#~ t,v=generateValidationGroups('all',0.6)
#~ writeTable(normalizeTable(gen_table('all',v)),T.dataRoot()+'/coiledcoil/all_types')


#~ writeValidationCandidates('homodimeric_parallel_can','homodimeric_parallel_val',v)
#~ writeValidationCandidates('homodimeric_antiparallel_can','homodimeric_antiparallel_val',v)
#~ writeValidationCandidates('heterodimeric_parallel_can','heterodimeric_parallel_val',v)
#~ writeValidationCandidates('heterodimeric_antiparallel_can','heterodimeric_antiparallel_val',v)
#~ writeValidationCandidates('all','all_val',v)



#~ dataFileCreation('homodimeric_parallel_dat','homodimeric_parallel_val',target_seq = "LLLLLLLLLLLL",target_type = ("homo","parallel",2), method = "All", use_big = 0 )
#~ study = CCStudy(data = 'homodimeric_parallel_dat')
#~ s1,a = study.doStudy(True)


#~ dataFileCreation('homodimeric_antiparallel_dat','homodimeric_antiparallel_val',target_seq = "LLLLLLLLLLLL",target_type = ("homo","antiparallel",2), method = "All", use_big = 0 )
#~ study = CCStudy(data = 'homodimeric_antiparallel_dat')
#~ s2,a = study.doStudy(True)
  

#~ dataFileCreation('heterodimeric_parallel_dat','heterodimeric_parallel_val',target_seq = "LLLLLLLLLLLL",target_type = ("hetero","parallel",2), method = "All", use_big = 0 )
#~ study = CCStudy(data = 'heterodimeric_parallel_dat')
#~ s3,a = study.doStudy(True)
   

#~ dataFileCreation('heterodimeric_antiparallel_dat','heterodimeric_antiparallel_val',target_seq = "LLLLLLLLLLLL",target_type = ("hetero","antiparallel",2), method = "All", use_big = 0 )
#~ study = CCStudy(data = 'heterodimeric_antiparallel_dat')
#~ s4,a = study.doStudy(True)

#~ ## ALL
#~ dataFileCreation('all_dat','all_val',target_seq = "LLLLLLLLLLLL",target_type = ("homo","parallel",2), method = "All", use_big = 0 )
#~ study = CCStudy(data = 'all_dat')
#~ s5,a = study.doStudy(True)

#~ file_stats('all')
#~ file_stats('all_val')
#~ file_stats('homodimeric_parallel_can')
#~ file_stats('homodimeric_parallel_val')
#~ file_stats('homodimeric_antiparallel_val')
#~ file_stats('heterodimeric_parallel_val')
#~ file_stats('heterodimeric_antiparallel_val')

#~ print "Homo par", s1 
#~ print "Homo antipar",s2 
#~ print "Hetero par",s3 
#~ print "Hetero antipar",s4 
#~ print "All",s5 
    
