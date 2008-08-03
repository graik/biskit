#!/usr/bin/env csh

set NB_NODES=%(nb_nodes)i
${LAMHOME}/bin/mpirun -np ${NB_NODES} -O  -ssi rpi tcp \
${AMBERHOME}/exe/sander -O -i sim.inp -o sim.out -p top.parm -c sim_2a.rst -ref sim_2a.rst -r sim.rst -x sim.crd -v sim.vel
