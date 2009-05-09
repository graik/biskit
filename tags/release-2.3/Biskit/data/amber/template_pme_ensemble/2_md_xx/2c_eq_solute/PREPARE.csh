#!/usr/bin/env csh
## this script is executed by amber_ensembleMD.py and then deleted.
##
## Links cannot be checked into CVS, create them on the fly

ln -s ../2b_eq_solvent_coP/sim.rst ./sim_2b.rst
ln -s ../../top.parm .
