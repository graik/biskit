from Biskit.polysys.coiledcoil import CoiledCoil
import sys

if __name__ == '__main__':
    
    ## Check if arguments are ok
    if len(sys.argv)<3:
        print "\nUsage: exe arg1 arg2\n"
        print "arg1->  File containing data to compare."
        print "arg2->  File containing scores."
        exit( -1 )
    
    ## Open data files...
    try:
        ccdata = open(sys.argv[1],"r")
        ccdata.close()
        ccdata = open(sys.argv[2],"r")
        ccdata.close()
    except (IOError,EOFError):
        print "-Error opening "+sys.argv[1]+": file doesn't exist."
        exit( -1 )
    
        
    c = CoiledCoil(sys.argv[2])
    
    c.parseData(sys.argv[1])
    
    ## Statistics
    sources = ["CC","Pa","Pair","Paper","So"]
    hits = {}
    tries = {}
    fails = {}
    for k in sources:
        hits[k] = 0
        tries[k] = 0
    sources.append("seq")
    
    
    heptads = {}
    for d in c.data.keys():
        ## Extract heptads
        heptads[d]={}
        for k in c.data[d].keys():
            if k in sources:
                heptads[d][k] = c.data[d][k]
    
    
    
    for h in heptads.keys():
        score = {}
        #~ print h
        there = heptads[h].keys()
        there.remove("seq")
        for k in there:
            com_keys = heptads[h].keys()
            com_keys.remove(k)
            com_keys.remove("seq")
            score[k] = 0
            for k2 in com_keys:
                if k != "Paper":
                    one = heptads[h][k]
                else:
                    one = heptads[h][k][0]
                if k2 != "Paper":
                    the_other = heptads[h][k2]
                else:
                    the_other = heptads[h][k2][0]
                
                if c.sameRegister(one,the_other,heptads[h]["seq"]):
                    score[k] += 1
        
        #~ print score    
        priorities = {"CC":4,"Pa":2,"Pair":4,"Paper":5,"So":3}

        best = (score.keys()[0],0)
        for k2 in score.keys():
            if best[1] < score[k2]:
                best = (k2,score[k2])
            if best[1] == score[k2]: ## In this case we have to manage priorities
                if priorities[k2] > priorities[best[0]]:
                    best  = (k2,score[k2])
        
        heptads[h]["best"] = best
        #~ print best
    
    for h in heptads.keys():
        there = heptads[h].keys()
        there.remove("seq")
        there.remove("best")
        for k in there:
            tries[k] +=1
            
            one = heptads[h][k]
            if k == "Paper":
                one = heptads[h][k][0]
                
            the_other = heptads[h][heptads[h]["best"][0]]
            if heptads[h]["best"][0] == "Paper":
                the_other = heptads[h][heptads[h]["best"][0]][0]
            
            if c.sameRegister(one,the_other,heptads[h]["seq"]):
                hits[k] += 1
        
    scores = {}
    for k in hits.keys():
        scores[k] = float(hits[k])/tries[k]
    
    
    print scores
    
    
    ## Get Alignment values
    
    
    target_hep = heptads["TARGET"][heptads["TARGET"]["best"][0]]
    if heptads["TARGET"]["best"][0]== "Paper":
        target_hep = heptads["TARGET"][heptads["TARGET"]["best"][0]][0]
    target = heptads["TARGET"]["seq"]
    target_reg=c.getRegister(target_hep,target)
    
    
    cross_scores = {}
    
    for h in heptads.keys():
        print h
        other_hep = heptads[h][heptads[h]["best"][0]]
        if heptads[h]["best"][0]== "Paper":
            other_hep = heptads[h][heptads[h]["best"][0]][0]
        other = heptads[h]["seq"]
        other_reg=c.getRegister(other_hep,other)
        print "l",len(other_reg),len(other)
        #~ print other, target, other_reg, target_reg
        
        res = c.tryAlignment(other, target, other_reg, target_reg)
        
        print res
        
        cross_scores[h]={}
        for i in ["heptad","res_like","charges"]:
            for j in ["heptad","res_like","charges"]:
                print res,res[j],res[j][1]
                cross_scores[h][(i,j)]= cross_scores[h][(j,i)] = c.getFitnessScore(j,other,target,res[j][1], other_reg,target_reg)
            
        print cross_scores[h]