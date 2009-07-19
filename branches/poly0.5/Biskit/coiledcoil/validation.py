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
            if contents[1] == '6':
                print l 
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
    
    file_par_hom = open("homodimeric_parallel","w")
    file_anti_hom = open("homodimeric_antiparallel","w")
    file_par_het = open("heterodimeric_parallel","w")
    file_anti_het = open("heterodimeric_antiparallel","w")
    
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


from coiledutils import getRegister
def gen_table(path , not_allowed=[]):
    file = open(path)
    lineas = file.readlines()
    file.close()
    lineas = [ l.strip() for l in lineas ]
    
    table ={}

    for l in lineas[1:]:
        contents = l.split()
        if len(contents) > 4 and len(contents[7])==7:
            register = getRegister(contents[7],contents[6])
            for i in range(len(contents[7])):  
                c = contents[7][i]
                r = register[i]
                if c in table:
                    table[c][r] += 1
                else:
                    table[c] = {'a':0,'b':0,'c':0,'d':0,'e':0,'f':0,'g':0}
                    table[c][r] = 1
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
    

file_stats('coils.db')   
splitInFiles('coils.db') 
file_stats('homodimeric_parallel')
print normalizeTable(gen_table('homodimeric_parallel'))
writeTable(normalizeTable(gen_table('homodimeric_parallel')),'table')

file_stats('homodimeric_antiparallel')
file_stats('heterodimeric_parallel')
file_stats('heterodimeric_antiparallel')


#~ def gen_table(file):
    
    
    
