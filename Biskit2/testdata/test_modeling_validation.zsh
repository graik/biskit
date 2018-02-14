#!/usr/bin/env zsh

echo " ***************************************************"
echo " ****** EXAMPLE SCRIPT FOR MODEL VALIDATION ********"
echo " ******       execute in ~biskit/test       ********"
echo " ***************************************************\n"

echo "In addition to the minimal requirements of Biskit this"
echo "script also utilize the following applocations that"
echo "need to be installed and in your executable path:"
echo " - ncbi-tools"
echo " - T-Coffee"
echo " - Modeller\n"

echo "In addition you will need to have pvm installed and the deamon running."


echo "RUNNING CALCULATIONS (will take approximately 20 minutes)\n"

cd Mod/project

echo " ============= Setup validation folder structure ============= " 
setup_validation.py


echo " ============= Align ============= " 
align_parallel.py


echo " ============= Modell ============= " 
model_parallel.py


echo " ============= Calculate benchmark data ============= "
benchmark.py

echo " ============= Analysis ============= "
analyze.py



echo "DONE\n"
