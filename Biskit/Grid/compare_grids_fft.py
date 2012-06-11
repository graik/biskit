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
#A = G1.trim(G2)
#A.source = G1.source
#B = G2.trim(G1)
#B.source = G2.source
A=G1
B=G2

#A.expand(6.)
#B.expand(6.)

mpos=A.data>0.0
A.data[mpos]=0.0
mpos=A.data<0.0
A.data[mpos]=1.0

mpos=B.data>0.0
B.data[mpos]=0.0
mpos=B.data<0.0
B.data[mpos]=1.0

#calculate reference max value

Fa=N.fft.rfftn(A.data)
Fb=N.fft.rfftn(A.data)

Fc=Fa.conj()*Fb

correl=N.fft.irfftn(Fc)

ref=correl.max()/N.sum(A.data>0.0)

#now compare:

Fa=N.fft.rfftn(A.data)
Fb=N.fft.rfftn(B.data)

Fc=Fa.conj()*Fb

correl=N.fft.irfftn(Fc)

npts=N.array([N.sum(A.data>0.1),N.sum(B.data>0.1)])

similarity = correl.max()/N.max(npts)

print "similarity : %.3f"%(similarity)