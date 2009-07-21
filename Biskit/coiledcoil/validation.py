import Biskit.tools as T
from coiledutils import getRegister
import random
from getcoil import dataFileCreation
from choosecoil import CCStudy
from Biskit.molUtils import allAA

## Statistics
def file_stats(path):
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

def splitInFiles(path):
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
        if len(contents) > 4 and contents[1] =='2':
            
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
    file = open(path)
    lineas = file.readlines()
    file.close()
    lineas = [ l.strip() for l in lineas ]
    
    table ={}
    for aa in allAA():
        table[aa] = {'a':0,'b':0,'c':0,'d':0,'e':0,'f':0,'g':0}
    
    print table
    for l in lineas[1:]:
        contents = l.split()
        if len(contents) > 4 and len(contents[7])==7 and not contents[0] in exclusion_list:
            register = getRegister(contents[7],contents[6])
            for i in range(len(contents[7])):  
                c = contents[7][i]
                r = register[i]
                table[c][r] += 1
                
    return table

def gen_table_samecontrib(path , exclusion_list=[]):
    """
    Same contribution.
    """
    
    file = open(path)
    lineas = file.readlines()
    file.close()
    lineas = [ l.strip() for l in lineas ]
    
    table ={}
    for aa in allAA():
        table[aa] = {'a':0,'b':0,'c':0,'d':0,'e':0,'f':0,'g':0}
    
    print table
    for l in lineas[1:]:
        contents = l.split()
        if len(contents) > 4 and len(contents[7])==7 and not contents[0] in exclusion_list:
            register = getRegister(contents[7],contents[6])
            for i in range(len(contents[7])):  
                c = contents[7][i]
                r = register[i]
                table[c][r] += 1./len(contents[7])
            
    return table

def gen_table_samecontrib_no_redundant(path , exclusion_list=[]):
    """
    Same contribution.
    """
    
    file = open(path)
    lineas = file.readlines()
    file.close()
    lineas = [ l.strip() for l in lineas ]
    
    table ={}
    for aa in allAA():
        table[aa] = {'a':0,'b':0,'c':0,'d':0,'e':0,'f':0,'g':0}
    
    print table
    for l in lineas[1:]:
        contents = l.split()
        if len(contents) > 4 and len(contents[7])==7 and not contents[0] in exclusion_list:
            register = getRegister(contents[7],contents[6])
            for i in range(len(contents[7])):  
                c = contents[7][i]
                r = register[i]
                if contents[8] == 'hetero':
                    table[c][r] += 1./len(contents[7])
                else:
                    table[c][r] += 1./(len(contents[7])*int( contents[1]))
            
    return table

def normalizeTable(table):
    reg = 'abcdefg'
    for r in reg:
        total = 0
        for aa in table:
            total += table[aa][r]
        
        for aa in table:
            table[aa][r] = table[aa][r]/float(total)

    return table        
    
    
def writeTable(table,path):
    file = open(path,"w")
    for aa in table:
        line = aa
        
        for r in table[aa]:
            line += " %.4f"%table[aa][r]
        
        file.write(line+"\n")
    file.close()


def generateValidationGroups(path,train_percent):
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
    file2 = open(path2,"w")
    
    file1 = open(path1)
    lineas = file1.readlines()
    file1.close()
    lineas = [ l.strip() for l in lineas ]
    file2.write(lineas[0]+"\n")
    pdbs = []
    for l in lineas[1:]:
        contents = l.split()
        if len(contents) > 4 and len(contents[6])>=length:
             file2.write(l+"\n")
    file2.close()

    
def writeValidationCandidates(path1,path2,val=[]):
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
    

import os
## All files
file_stats('coils.db')   


lengthFilter('coils.db','coils_filtered.db',0)
#~ file_stats('coils_long.db')  
#~ splitInFiles('coils_long.db') 

os.system("cp coils_filtered.db all")
splitInFiles('all') 


## homodimeric_parallel
file_stats('homodimeric_parallel_can')
t,v=generateValidationGroups('homodimeric_parallel_can',0.8)
writeTable(normalizeTable(gen_table_samecontrib_no_redundant('homodimeric_parallel_can',v)),T.dataRoot()+'/coiledcoil/homodimeric_parallel')
writeValidationCandidates('homodimeric_parallel_can','homodimeric_parallel_val',v)
dataFileCreation('homodimeric_parallel_dat','homodimeric_parallel_val',target_seq = "LLLLLLLLLLLL",target_type = ("homo","parallel",2), method = "All", use_big = 0 )
study = CCStudy(data = 'homodimeric_parallel_dat')
s1,a = study.doStudy(True)
  

#~ file_stats('homodimeric_antiparallel_can')
#~ t,v=generateValidationGroups('homodimeric_antiparallel_can',0.8)
#~ writeTable(normalizeTable(gen_table('homodimeric_antiparallel_can',v)),T.dataRoot()+'/coiledcoil/homodimeric_antiparallel')
#~ writeValidationCandidates('homodimeric_antiparallel_can','homodimeric_antiparallel_val',v)
#~ dataFileCreation('homodimeric_antiparallel_dat','homodimeric_antiparallel_val',target_seq = "LLLLLLLLLLLL",target_type = ("homo","antiparallel",2), method = "All", use_big = 0 )
#~ study = CCStudy(data = 'homodimeric_antiparallel_dat')
#~ s2,a = study.doStudy(True)
  

#~ file_stats('heterodimeric_parallel_can')
#~ t,v=generateValidationGroups('heterodimeric_parallel_can',0.8)
#~ writeTable(normalizeTable(gen_table('heterodimeric_parallel_can',v)),T.dataRoot()+'/coiledcoil/heterodimeric_parallel')
#~ writeValidationCandidates('heterodimeric_parallel_can','heterodimeric_parallel_val',v)
#~ dataFileCreation('heterodimeric_parallel_dat','heterodimeric_parallel_val',target_seq = "LLLLLLLLLLLL",target_type = ("hetero","parallel",2), method = "All", use_big = 0 )
#~ study = CCStudy(data = 'heterodimeric_parallel_dat')
#~ s3,a = study.doStudy(True)
   

#~ file_stats('heterodimeric_antiparallel_can')
#~ t,v=generateValidationGroups('heterodimeric_antiparallel_can',0.8)
#~ writeTable(normalizeTable(gen_table('heterodimeric_antiparallel_can',v)),T.dataRoot()+'/coiledcoil/heterodimeric_antiparallel')
#~ writeValidationCandidates('heterodimeric_antiparallel_can','heterodimeric_antiparallel_val',v)
#~ dataFileCreation('heterodimeric_antiparallel_dat','heterodimeric_antiparallel_val',target_seq = "LLLLLLLLLLLL",target_type = ("hetero","antiparallel",2), method = "All", use_big = 0 )
#~ study = CCStudy(data = 'heterodimeric_antiparallel_dat')
#~ s4,a = study.doStudy(True)



#~ ## ALL
#~ print "ALL","\n\n\n\n\n\n\n\n\n\n\n"
#~ t,v=generateValidationGroups('all',0.8)
#~ writeTable(normalizeTable(gen_table('all',v)),T.dataRoot()+'/coiledcoil/all_types')
#~ writeValidationCandidates('all','all_val',v)
#~ dataFileCreation('all_dat','all_val',target_seq = "LLLLLLLLLLLL",target_type = ("homo","parallel",2), method = "All", use_big = 0 )
#~ study = CCStudy(data = 'all_dat')
#~ s5,a = study.doStudy(True)

file_stats('homodimeric_parallel_val')
#~ file_stats('homodimeric_antiparallel_can')
#~ file_stats('heterodimeric_parallel_can')
#~ file_stats('heterodimeric_antiparallel_can')
#~ file_stats('all')

print "Homo par", s1 
#~ print "Homo antipar",s2 
#~ print "Hetero par",s3 
#~ print "Hetero antipar",s4 
#~ print "All",s5 
    
    
