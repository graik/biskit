#!/usr/bin/env csh

set NB_NODES=%(nb_nodes)i

${LAMHOME}/bin/mpirun -np ${NB_NODES} -O -ssi rpi tcp \
${AMBERHOME}/exe/sander -O -i sim_0.inp -o sim_0.out -p top.parm -c sim_2b.rst -ref sim_2b.rst -r sim_0.rst -x sim_0.crd -v sim_0.vel

${LAMHOME}/bin/mpirun -np ${NB_NODES} -O -ssi rpi tcp \
${AMBERHOME}/exe/sander -O -i sim_1.inp -o sim_1.out -p top.parm -c sim_0.rst -ref sim_0.rst -r sim_1.rst -x sim_1.crd -v sim_1.vel

${LAMHOME}/bin/mpirun -np ${NB_NODES} -O -ssi rpi tcp \
${AMBERHOME}/exe/sander -O -i sim_2.inp -o sim_2.out -p top.parm -c sim_1.rst -ref sim_1.rst -r sim_2.rst -x sim_2.crd -v sim_2.vel

${LAMHOME}/bin/mpirun -np ${NB_NODES} -O -ssi rpi tcp \
${AMBERHOME}/exe/sander -O -i sim_3.inp -o sim_0.out -p top.parm -c sim_2.rst -ref sim_2.rst -r sim_3.rst -x sim_3.crd -v sim_3.vel

${LAMHOME}/bin/mpirun -np ${NB_NODES} -O -ssi rpi tcp \
${AMBERHOME}/exe/sander -O -i sim_4.inp -o sim_4.out -p top.parm -c sim_3.rst -ref sim_3.rst -r sim_4.rst -x sim_4.crd -v sim_4.vel
