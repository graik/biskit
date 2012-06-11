#!/usr/bin/python
import Biskit.Grid.gridManager as gr
import numpy as N
import sys
import pylab

__doc__ = """USE: python compareAllTanimoto.py out_file_name referencePDB [grid list]"""

if len(sys.argv)< 4:
  sys.exit("USE: python compareAllTanimoto.py out_file_name referencePDB [grid list]")


outFile = sys.argv[1]
referencePDB = sys.argv[2]
fileList= sys.argv[3:]
fileList = N.array(fileList)

gridList = [ gr.get_grid(file) for file in fileList ]

[grid.update(grid.mergeConserveProt(referencePDB, 3.5)) for grid in gridList]

result = N.zeros([fileList.shape[0], fileList.shape[0] ])

for id1, grid1 in enumerate(gridList):
    if id1 == len(gridList)-1 :break
    for id2, grid2 in enumerate(gridList[id1+1:]):
        result[id1,id2+id1+1] = gr.tanimotoIdentityCutoff(grid1, grid2, cutoff=10.)

resultWTitle = N.vstack((fileList, result))

#result = result + result.T + (N.eye(result.shape[0])*(result.max()+0.05))
result = result + result.T + N.eye(result.shape[0])

pylab.matshow(result)
pylab.title(outFile)
pylab.clim(vmin=0,vmax=1)
pylab.colorbar()
pylab.savefig(outFile+'.png')
#pylab.show()
N.savetxt(outFile+'.txt',resultWTitle,fmt='%s', delimiter=' ')
