from residuepicker import residuePicker
from folderdatacreator import FoldDataCreator
import os
from restools import expandedPoints,linearize
import numpy as N
from vectors import normalized,sub
from quaternion import rotquat, rotmat
from Biskit import PDBModel
from vectors import matrix2list
import random

FILL_DATABASES = False
MAX_ITERATIONS = 1500
MAX_FILES = 50


## Create a new database
rpicker = residuePicker()
## Create a new folder tool
fold = FoldDataCreator()

if FILL_DATABASES:
    ## Which pdbs?
    files = os.listdir("./data/pdb")
    ## Choose MAX_FILES randomly
    rand = random.sample(files,MAX_FILES)
    ## Fill the DB
    fails = 0
    for f in range(MAX_FILES):
        print "Extracting...",rand[f]
        try:
            rpicker.extractFromPDB('./data/pdb/'+rand[f])
        except:
            fails = fails+1
    rpicker.save()
    
    rand = random.sample(files,MAX_FILES)
    for f in range(MAX_FILES):
        print "Processing...",rand[f]
        try:
            fold.process_protein('./data/pdb/'+rand[f])
            fold.consumeData()
        except:
            fails = fails+1
    fold.save()


    print "Fails:",fails
else:
   ## Start folding of some structures
    _1AFP_sequence = "ATYNGKCYKKDNICKYKAQSGKTAICKCYVKKCPRDGAKCEFDSYKGKCYC"
    sequences = [_1AFP_sequence]
    
    
    for s in sequences:
        points = expandedPoints(len(s))
        chain = rpicker.createChain(s,points)
        chain.writePdb('aux.pdb')
        chain = PDBModel('aux.pdb')
        chain_ant = chain.clone()
        
        for i in range(MAX_ITERATIONS):
            
            ## Predict next movements
            translations = fold.getTranslations(chain)
            
            print "Applying new translation vectors..."
            print len(translations)
            print translations
            ##Calculate newpoints based on translation and rotation
            for j in range(len(s)-1):
                print s[j]
                ## Center the points in the one of interest
                print points[j]
                aux = N.array([points[j][0],points[j][1],points[j][2]])
                for k in range(len(points)):
                    points[k] = points[k] - aux

                ## Do rotation
                v1 = points[j+1]
                v2 = [1,0,0]
                
                R= linearize(None,v1,v2)
                
                for k in range(len(points)):
                    points[k] = N.array(matrix2list(points[k] * R))
                
                ## Do the translation
                points[j+1]=[translations[j][0],translations[j][1],translations[j][2]]
                
                for k in range(j+2,len(points)):
                    points[k] = points[k] + sub(points[j+1],translations[j])
                
            ## Center to the first and rotate to the base position
            
            aux = N.array([points[0][0],points[0][1],points[0][2]])
            for k in range(len(points)):
                points[k] = points[k] - aux
            v1 = points[1]
            v2 = [1,0,0]
            v1a = N.array([normalized(v1),[0,0,1],[0,0,1]])
            v2a = N.array([normalized(v2),[0,0,1],[0,0,1]])
            q  = rotquat(v1a,v2a)
            R = N.transpose(N.matrix(rotmat(q[0])))
            print R
            for k in range(len(points)):
                points[k] = N.array(matrix2list(points[k] * R))
                    
            ## Create next conformation
            chain = rpicker.createChain(s,points)
            chain.writePdb('aux.pdb')
            chain = PDBModel('aux.pdb')
            if i % 5 == 0:
                chain.writePdb("./results/frame"+str(i)+".pdb")
            if chain.rms(chain_ant)< 0.01:
                i = MAX_ITERATIONS
            chain_ant = chain.clone()
