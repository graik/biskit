#!/usr/bin/python2.6
import GRID as G
import numpy as N
import sys, os

#if len(sys.argv)<4:
  #sys.exit("USAGE: compare_grids.py GRID_file1 GRID_file2 Output_File_name")

INVERSE = 0
EXPAND = 1

#HERE FILE NAMES OF THE TWO GRIDS TO COMPARE (should be Free Energy grids) and the output
grid1_file = sys.argv[1]
grid2_file = sys.argv[2]
CUTPLIST = [60., 40., 20., 10.]
#out = sys.argv[3]   #Output log name

print "\n=================================================================="

#----------------------------- RAW GRID comparison first -------------------------------------
print "Loading grids and trimming...",
#Load Grids
G1 = G.get_grid(grid1_file)
G2 = G.get_grid(grid2_file)

#Trim grids to obtain common space
A = G1.trim(G2)
A.source = G1.source
B = G2.trim(G1)
B.source = G2.source
print "DONE"
print "GRID A: ",A.source
print "GRID B: ",B.source
#AB_out = str(os.path.splitext(out)[0])+"_ALL"+str(os.path.splitext(out)[1])

# ----------------- TANIMOTO COEFFICIENT CALCULATION FOR RAW GRIDS --------------------------------
#Convert grids to bits giving a cutoff of 40% the most negative points (60% minimum value)

#CUTOFF = (100. - CUTP) / 100.

#A_cut = A.data.min()*CUTOFF
#B_cut = B.data.min()*CUTOFF

for CUTP in CUTPLIST:
    #Flatten negative data, sort it and get the value of the 20% position of negative point as cutoff
    fsort_A = N.array( A.data[A.data<0].flat )
    fsort_A.sort()
    CUT_INDX = int( fsort_A.shape[0] * (CUTP/100.) )
    A_CUT_VAL = fsort_A[CUT_INDX]
    
    fsort_B = N.array( B.data[B.data<0].flat )
    fsort_B.sort()
    CUT_INDX = int( fsort_B.shape[0] * (CUTP/100.) )
    B_CUT_VAL = fsort_B[CUT_INDX]
    
    #Mask TRUE (1) values under cutoff or between cutoff and zero if INVERSE = 1
    if INVERSE == 1:
        A_m1 = (A.data > A_CUT_VAL) * (A.data < 0)
        B_m1 = (B.data > B_CUT_VAL) * (B.data < 0)
    else:
        A_m1 = A.data <= A_CUT_VAL  
        B_m1 = B.data <= B_CUT_VAL
    
    
    #EXPAND POSITIVE POINTS 0.5 ANGSTROMS AROUND
    if EXPAND == 1:
        E_A = A.copy()
        E_B = B.copy()
        E_A.data = A_m1*1
        E_A.expand(2)   #Expand to avoid bad indexing
        E_B.data = B_m1*1
        E_B.expand(2)
        for grid in (E_A, E_B):
            point_lst = N.vstack(N.nonzero(grid.data==1)).T
            for point in point_lst:
                grid.cancelpoints(point = point, cutoff = 0.5, value = 1)
        E_A.contract(2)
        E_B.contract(2)
        A_m1 = E_A.data == 1
        B_m1 = E_B.data == 1
    
    #Mix bitwise AND
    MM = A_m1 * B_m1
    
    #Tanimoto calculation
    Na, Nb, Nab = len(N.nonzero(A_m1)[0]), len(N.nonzero(B_m1)[0]), len(N.nonzero(MM)[0])
    Tc = float(Nab) / (Na + Nb - Nab)
    
    print "\n-------------------------------------------------"
    print "CUT PERCENTAGE: ", CUTP
    print "Grid A MIN and CUTOFF: ",A.data.min(), A_CUT_VAL
    print "Grid B MIN and CUTOFF: ",B.data.min(), B_CUT_VAL
    print "Na:",Na,"Nb:",Nb,"Nab",Nab
    print "TANIMOTO COEFFICIENT: ", Tc
    
print "\n==================================================================="
# ----------------- HIERARCHICAL CLUSTERING ------------------------------
import scipy.cluster.hierarchy as HC

for (grid,cut) in [(A, A_CUT_VAL),(B, B_CUT_VAL)]:
    GOOD = N.vstack(N.where(grid.data <= cut)).T
    for i, el in enumerate(GOOD):
        GOOD[i] = grid.getCartesian(el)
    
    LINK = HC.linkage(GOOD, method='average')
    
    HC.matplotlib.pylab.figure()
    DEND = HC.dendrogram(LINK)
    
    for t in (8, 10, 12):
        T = HC.fcluster(LINK, t, criterion='distance')
        f = open('%s_clusterpdb_%i.pdb'%(grid.source,t),'w')
        for i,pos in enumerate(GOOD):
            f.write('ATOM  %5d%4s  PTH  %4d   %8.3f%8.3f%8.3f%6.2f%6.2f\n'%(i,"C",T[i],pos[0],pos[1],pos[2],0.0,0.0))
        f.close()

HC.matplotlib.pylab.show()
## ----------------- Generate Grids with the negative points only ----------------------------------
##Masking negative points
#m_base = A.data<0
#m_check = B.data<0

##Merge both masks 
#MASK = m_base * m_check

#mA = A.copy()
#mA.update(mA.data[MASK])
#mB = B.copy()
#mB.update(mB.data[MASK])
#mAmB_out = str(os.path.splitext(out)[0])+"_NEG"+str(os.path.splitext(out)[1])

## ----- We will reescale all grid points with mean values of the positive and negative side as 1, -1 -----------

#def scale(val, negative_mean):
  #"""Scale negative values with negative_mean"""
  #if val >0:
    #return val
  #else: 
    #return val/negative_mean

#def negmean(grid):
  #"""Return mean value of the negative points"""
  #neg_m = grid.data<0
  #negative_mean = grid.data[neg_m].mean()
  #return -1*negative_mean

#scaling = N.vectorize(scale)

#SA = A.copy()

#m1_A = SA.data>1    #Replace all numbers >1 with 1
#SA.data[m1_A] = 1
#neg_mean = negmean(SA)
#SA.update( scaling(SA.data, neg_mean) )

#SB = B.copy()

#m1_B = SB.data>1    #Replace all numbers >1 with 1
#SB.data[m1_B] = 1
#neg_mean = negmean(SB)
#SB.update( scaling(SB.data, neg_mean) )
#SASB_out = str(os.path.splitext(out)[0])+"_SCALED"+str(os.path.splitext(out)[1])

##Transform all points >1 to 1


#print "BASE: ",G1.source
#print "CHECK: ",G2.source

##-----------------     Statistical analysis and file writing for all grids      -----------------------------------

#for (base, check, fout, desc) in [(A, B, AB_out, "ALL"), (mA, mB, mAmB_out, "NEGATIVE"), (SA, SB, SASB_out, "SCALED")]:
    
    #print "Calculating %s points stats..."%(desc)
    
    #Es = base.data - check.data         #Error array with sign
    #Ens = abs(base.data - check.data)   #Error array without sign
    
    ##Number of points with error >0.5
    #epoints = len(N.where(Ens>0.5)[0])
    #r_coef = N.corrcoef(base.data.flat, check.data.flat)[0,1]
    
    #fo = open(fout, 'w')
    
    #fo.write("\t\t## GRIDs COMPARISON ##\n\nBase: %s\tCheck: %s\n\n--------------- Statistics %s POINTS -----------------\n\n\tError with sign :\tError without sign:\n\nMeans:\t\t%9.5f\t%9.5f\nStd Dev.:\t%9.5f\t%9.5f\nMin.:\t\t%9.5f\t%9.5f\nMax.:\t\t%9.5f\t%9.5f\n\nPercentage of points with Error >0.5: %5.2f \n\nSquared Correlation Coeficient: %7.4f\n________________________________________________________\n\n"%(G1.source, G2.source, desc, Es.mean(), Ens.mean(), Es.std(), Ens.std(), Es.min(), Ens.min(), Es.max(), Ens.max(), epoints/float(Ens.size)*100, r_coef**2))
    
    #print "Writing data..."
    #for x in xrange(base.data.size):
        #fo.write("%8.3f\t%8.3f\t%8.3f\t%8.3f\n"%(base.data.flat[x],check.data.flat[x], Es.flat[x], Ens.flat[x]))
    
    #fo.close()
    #print "DONE"
