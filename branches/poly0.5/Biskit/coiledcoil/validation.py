## Statistics
def file_stats(path):
    file = open(path)
    lineas = file.readlines()
    file.close()
    lineas = [ l.strip() for l in lineas ]
    
    pdbs = {}
    coil_oligo = {'2':0,'3':0,'4':0,'5':0,'6':0,'7':0,'8':0}
    paralel_homo = 0
    paralel_hetero = 0
    antip_homo = 0 
    antip_hetero = 0
    seq_length = 0
    
    for l in lineas:
        contents = l.split()
        if contents[0] in pdbs:
            pdbs[contents[0]] += 1./int(contents[1])
        else:    
            pdbs[contents[0]] = 1./int(contents[1])
        
        coil_oligo[contents[1]] += 1./int(contents[1])
        
        seq_length += (1./int(contents[1])) * len(contents[6])
        
        if contents[2] == 'parallel' and contents[8] == 'homo' :
            paralel_homo += 1./int(contents[1])
        elif contents[2] == 'parallel' and contents[8] == 'hetero':
            paralel_hetero += 1./int(contents[1])
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
        
    print "Total Coils (1):",total_coils
    print "Total Coils (2):",total_coils2
    print "Total Files:",total_files
    print "Coils per file:",total_coils / float(total_files)
    print "Oligomerization:",coil_oligo
    print "Mean sequence length:", seq_length / total_coils
    
    
    
#~ def gen_table(file):
    
    
    