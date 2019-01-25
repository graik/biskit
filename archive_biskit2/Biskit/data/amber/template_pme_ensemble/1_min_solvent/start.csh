#!/usr/bin/env csh

set NB_NODES=%(nb_nodes)i
${LAMHOME}/bin/mpirun -np ${NB_NODES} -O -ssi rpi tcp \
${AMBERHOME}/exe/sander -O -i min.inp -o min.out -p top.parm -c 0.crd -ref 0.crd -x min.crd -r min.rst
