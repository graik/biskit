#!/usr/bin/env zsh

echo " ****************************************************"
echo " ******** EXAMPLE SCRIPT FOR A DOCKING RUN **********"
echo " ********     execute in ~biskit/test      **********"
echo " ****************************************************\n"

echo "In addition to the minimal requirements of Biskit this"
echo "script also utilize the following applocations that"
echo "need to be installed and in your executable path:"
echo " - ncbi-blast (http://www.ncbi.nlm.nih.gov/blast/)"
echo " - Xplor-nih (http://nmr.cit.nih.gov/xplor-nih/)"
echo " - HEX rigid body docking (http://www.csd.abdn.ac.uk/hex/)\n"

echo "Additional applications that will be used if they are available:"
echo " - WhatIf"
echo " - Hmmer"
echo " - SurfaceRacer"
echo " - Fold-X\n"

echo "RUNNING CALCULATIONS (will take approximately 40 minutes)\n"

echo " ============= Cleaning up receptor pdb file ============= " 
cd rec
pdb2xplor.py -i 1A2P_rec_original.pdb 
xplor-nih < 1A2P_rec_original_generate.inp > 1A2P_rec_original_generate.log
echo "DONE cleaning\n"

echo " ============= Cleaning up ligand pdb file ============= " 
cd ../lig
pdb2xplor.py -i 1A19_lig_original.pdb -c 1
xplor-nih < 1A19_lig_original_generate.inp > 1A19_lig_original_generate.log
echo "DONE cleaning\n"

echo " ============= Cleaning up reference complex pdb file ============= " 
cd ../com
pdb2xplor.py -i 1BGS_edited.pdb
xplor-nih < 1BGS_edited_generate.inp  > 1BGS_edited_generate.log
pdb2complex.py -c 1BGS.pdb  -r 0 -l 1
echo "DONE cleaning\n"

echo " ============= Creating docking project root ============= "
mkdir ../dock
cd ../dock

echo " ============= Writing receptor files to be used for docking ============= "
mkdir rec
cd rec
PCR2hex.py -psf ../../rec/1A2P.psf -pdb ../../rec/1A2P.pdb
dope.py -s ../../rec/1A2P.pdb -so ../../rec/1A2P_dry.model -i 1A2P.model -dic 1A2P_model.dic
echo "DONE\n"

echo " ============= Writing ligand files to be used for docking ============= "
mkdir ../lig
cd ../lig
PCR2hex.py -psf ../../lig/1A19.psf -pdb ../../lig/1A19.pdb
dope.py -s ../../lig/1A19.pdb  -so ../../lig/1A19_dry.model -i 1A19.model  -dic 1A19_model.dic
echo "DONE\n"

echo " ============= Writing reference complex to be used for docking ============= "
mkdir ../com
cd ../com
PCR2hex.py -psf ../../com/1BGS.psf -pdb ../../com/1BGS.pdb
dope.py -s ../../com/1BGS.pdb  -so ../../com/1BGS_dry.model -i 1BGS.model  -dic 1BGS_model.dic
echo "DONE\n"

echo " ============= Writing HEX docking macro file ============= "
mkdir ../hex
cd ../hex
hexInput.py -r ../rec/1A2P_hex.pdb -l ../lig/1A19_hex.pdb -c ../com/1BGS_hex.pdb
echo "DONE\n"

echo " ============= Running HEX docking ============= "
if [ `cat /proc/cpuinfo | grep processor | wc -l` -ge 2  ]; then
	echo "Running HEX docking on 2 CPUs"
	hex -ncpu 2 < 1A2P-1A19_hex.mac > 1A2P-1A19_hex.log
else
	echo "Running HEX docking on 1 CPU"
	hex < 1A2P-1A19_hex.mac > 1A2P-1A19_hex.log
fi
# wait until hex is done before proceeding
while [[ `ps aux | grep $HEX_ROOT/exe/hex | wc -l` -ge 2 ]]; do
	sleep 100
done
echo "DONE docking\n"

echo " ============= Parsing HEX output ============= "
hex2complex.py -rec ../rec/1A2P_model.dic -lig ../lig/1A19_model.dic -hex 1A2P-1A19_hex.out -p
echo "DONE\n"

echo " ============= Calculateing fractions of native contacts ============= "
contacter.py -i complexes.cl -ref ../../com/ref.complex -a
inspectComplexList.py complexes.cl
echo "DONE\n"
