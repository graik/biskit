#!/usr/bin/env csh

set NB_NODES=%(nb_nodes)i

${LAMHOME}/bin/mpirun -np ${NB_NODES} -O -ssi rpi tcp \
${AMBERHOME}/exe/sander -O -i sim.inp -o sim.out -p top.parm -c eq.rst -r sim.rst -x sim.crd -v sim.vel

echo "Post-processing..."

ptraj top.parm ../../postprocess.ptraj

gzip sim.crd
gzip sim.vel&

amber2traj.py -i sim_protein.crd -o traj_protein.dat -r ../../0_dry.pdb -b -wat -hyd -rnres

gzip sim_protein.crd

echo "Done"
