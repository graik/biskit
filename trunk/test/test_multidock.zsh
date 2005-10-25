#!/usr/bin/env zsh

echo " ****************************************************"
echo " ***** EXAMPLE SCRIPT FOR A MULTI-DOCKING RUN *******"
echo " *****   execute in ~biskit/test, requires    *******"
echo " *****    that test_docking has been run      *******"
echo " ****************************************************\n"

echo "In addition to the minimal requirements of Biskit this"
echo "script also utilize the following applocations that"
echo "need to be installed and in your executable path:"
echo " - Xplor capable of pcr restraints (avaliable upon request)"
echo " - HEX rigid body docking (http://www.csd.abdn.ac.uk/hex/)\n"

echo "Additional applications that will be used if they are available:"
echo " - WhatIf"
echo " - Hmmer"
echo " - SurfaceRacer"
echo " - Fold-X\n"

echo "RUNNING CALCULATIONS (will take approximately 60-80 minutes)\n"
echo "REMINDER: Check that pvm is running.\n"

date

echo " ============= Run a short ligand MD on the localhost ============= " 
# geting test root path
root=`echoTestRoot.py | egrep "^/.*"`

mkdir lig_pcr_00
runPcr.py -t lig -r lig_pcr_00 -h localhost -nstep 10
while [[ `ps aux | egrep  "$root/lig_pcr_00/pcr_00" | wc -l` -ge 2 ]]; do
	sleep 100
done
echo "MD DONE\n"

echo " ============= Convert snapshots to models ============= " 
cd lig_pcr_00
mkdir model
pdb2model.py -i pcr_00/*.pdb -a -wat -o model/
echo "DONE\n"

echo " ============= Create a trajectory from models ============= " 
pdb2traj.py -i model/*.model -r ../lig/1A19.pdb -wat -f
thinTraj.py -i traj.dat
echo "DONE\n"

cd ..
mkdir multidock
mkdir multidock/lig
cd multidock/lig

echo " ============= Cluster ligand trajectory ============= " 
selectModels.py -i ../../lig_pcr_00/traj.dat -n 2
echo "DONE clustering\n"

echo " ============= Prepare cluster centers for docking ============= " 
PCR2hex.py  -psf ../../lig/*psf -pdb *model
echo "DONE\n"

echo " ============= Add structure realted data to models ============= " 
dope.py -s ../../lig/1A19.pdb  -so ../../lig/1A19_dry.model -i *.model  -dic 1A19_model.dic
echo "DONE\n"
  
cd ..
ln -s ../dock/rec .
ln -s ../dock/com .

echo " ============= Start multi-docking ============= " 
multidock.py -rdic rec/1A2P_model.dic -ldic lig/1A19_models.dic -com com/1BGS_hex.pdb -hex
echo "DONE docking\n"

cd hex*
contacter.py -i complexes.cl -a -ref ../../com/ref.complex
