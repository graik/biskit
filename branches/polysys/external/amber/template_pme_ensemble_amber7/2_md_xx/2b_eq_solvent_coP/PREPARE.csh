#!/usr/bin/env csh
## this script is executed by amber_ensembleMD.py and then deleted.
##
## Links cannot be checked into CVS, create them on the fly

ln -s ../2a_eq_solvent_coV/sim.rst ./sim_2a.rst
ln -s ../../top.parm .
