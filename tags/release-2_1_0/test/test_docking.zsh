#!/usr/bin/env zsh

echo " ****************************************************"
echo " ******** EXAMPLE SCRIPT FOR A DOCKING RUN **********"
echo " ********     execute in ~biskit/test      **********"
echo " ****************************************************\n"

echo "In addition to Biskit this script also utilizez the following"
echo "applications that need to be installed and in your executable path:"
echo " - ncbi-blast (http://www.ncbi.nlm.nih.gov/blast/)"
echo " - Xplor-nih (http://nmr.cit.nih.gov/xplor-nih/)"
echo " - HEX rigid body docking (http://www.csd.abdn.ac.uk/hex/)"
echo " - surfrace\n"

echo "The last part of the script calls contacter.py and needs PVM"
echo "installed and configured." 

echo "RUNNING CALCULATIONS  will take approximately 40 minutes\n"
echo "Use option 'dry' to skip any calls to external applications:"
echo "   test_docking.zsh dry"
echo
echo "This will only re-create pickled objects, except that it still calls"
echo "the surfrace application (if available) to add surface profiles."
echo "The Biskit/Dock/ContactMaster.py test case depends on these profiles."

echo "Use option 'clean' to remove intermediate files that are not"
echo "needed by the Biskit test cases and are not part of the repository.\n"

##############
## Preparation
##############

if [[ $1 = dry || $2 = dry ]]; then
    export dry=true
else
    export dry=false
fi

if [[ $1 = clean || $2 = clean ]]; then
    export clean=true
else
    export clean=false
fi

if [[ $dry = true ]]; then
    echo " ======== DRY-RUN without calling external applications====== "
    echo
fi

root=`pwd`
PATH=$PATH:`pwd`/../scripts/Biskit:`pwd`/../scripts/Dock:`pwd`/../scripts/analysis

##################################
## Call scripts and generate tests
##################################

echo " ============= Cleaning up receptor pdb file ============= " 
cd $root/rec
pdb2xplor.py -i 1A2P_rec_original.pdb 
if [[ $dry = false ]]; then
    xplor-nih < 1A2P_rec_original_generate.inp > 1A2P_rec_original_generate.log
fi
echo "DONE cleaning\n"

echo " ============= Cleaning up ligand pdb file ============= " 
cd $root/lig
pdb2xplor.py -i 1A19_lig_original.pdb -c 1
if [[ $dry = false ]]; then
    xplor-nih < 1A19_lig_original_generate.inp > 1A19_lig_original_generate.log
fi
echo "DONE cleaning\n"

echo " ============= Cleaning up reference complex pdb file ============= " 
cd $root/com
pdb2xplor.py -i 1BGS_edited.pdb
if [[ $dry = false ]]; then
    xplor-nih < 1BGS_edited_generate.inp  > 1BGS_edited_generate.log
fi    
pdb2complex.py -c 1BGS.pdb  -r 0 -l 1
echo "DONE cleaning\n"

echo " ============= Creating docking project root ============= "
mkdir $root/dock
cd $root/dock
echo "Done\n"

echo " ============= Writing receptor files to be used for docking ============= "
mkdir $root/dock/rec
cd $root/dock/rec
PCR2hex.py -psf ../../rec/1A2P.psf -pdb ../../rec/1A2P.pdb
dope.py -s ../../rec/1A2P.pdb -so ../../rec/1A2P_dry.model -i 1A2P.model -dic 1A2P_model.dic -p dens surf
echo "DONE\n"

echo " ============= Writing ligand files to be used for docking ============= "
mkdir $root/dock/lig
cd $root/dock/lig
PCR2hex.py -psf ../../lig/1A19.psf -pdb ../../lig/1A19.pdb
dope.py -s ../../lig/1A19.pdb  -so ../../lig/1A19_dry.model -i 1A19.model  -dic 1A19_model.dic -p dens surf
echo "DONE\n"

echo " ============= Writing reference complex pdb file ============= "
mkdir $root/dock/com
cd ../com
PCR2hex.py -psf ../../com/1BGS.psf -pdb ../../com/1BGS.pdb
dope.py -s ../../com/1BGS.pdb  -so ../../com/1BGS_dry.model -i 1BGS.model  -dic 1BGS_model.dic -p dens surf
echo "DONE\n"

echo " ========= Creating reference complex to be used for docking ========= "
pdb2complex.py -c $root/com/1BGS.pdb  -r 0 -l 1
echo "DONE\n"

echo " ============= Writing HEX docking macro file [needs surfrace] ======= "
mkdir $root/dock/hex
cd $root/dock/hex
if [[ $dry = false ]]; then
    hexInput.py -r ../rec/1A2P_hex.pdb -l ../lig/1A19_hex.pdb -c ../com/1BGS_hex.pdb
    echo "DONE\n"
else
    echo "Skipped\n"
fi

echo " ============= Running HEX docking [needs hex] ======== "
cd $root/dock/hex
if [[ $dry = false ]]; then
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
    echo "DONE\n"
else
    echo "Skipped\n"
fi


echo " ============= Parsing HEX output ============= "
hex2complex.py -rec ../rec/1A2P_model.dic -lig ../lig/1A19_model.dic -hex 1A2P-1A19_hex.out -p
echo "DONE\n"

echo " ============= Calculating fractions of native contacts ============= "
contacter.py -i complexes.cl -ref ../../com/ref.complex -a
inspectComplexList.py complexes.cl
a_compare_rms_vs_fnc.py -i complexes_cont.cl
echo "DONE\n"

############################
## Remove what is not needed
############################

if [[ $clean = true ]]; then
    echo "============ Clean up: keep only files needed for tests ========"
    cd $root
    rm com/1BGA_seg.PDB
    rm com/1BGB_seg.PDB
    rm com/1BGS_dry.model
    rm com/1BGS_edited.log
    rm com/1BGS_edited_generate.inp
    rm com/1BGS_waters.pdb
    rm -r dock/com
    rm dock/hex/complexes.eps
    rm dock/lig/1A19_hex.pdb
    rm dock/rec/1A2P_hex.pdb
    rm lig/1A19_lig_original.log
    rm lig/1A19_lig_original_generate.inp
    rm lig/1A19_waters.pdb
    rm lig/1A1B_seg.PDB
    rm rec/1A2A_seg.PDB
    rm rec/1A2P_rec_original.log
    rm rec/1A2P_rec_original_generate.inp
    rm rec/1A2P_waters.pdb
fi