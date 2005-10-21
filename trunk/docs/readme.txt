Biskit
======

scripts and data for structure analysis, flexible docking,
and homology modelling
Authors: Johan & Raik

created: 22/04/02

Installation
============

- Biskit
  check out the project from CVS server (creates folder biskit)
	% cvs co biskit

- Python (version 2.4 recomended)
  install from rpms python and python-devel (source / include files) 
   or 
  from source
  	- download from http://www.python.org/
	- in: python-2.4.x:  ./configure
                             make
		             make install  (installs in /usr/local)

- Plotutils
  install from rpm
  or
  from source
        - download from http://www.gnu.org/software/plotutils/
        - in plotutils-2.4.1:  ./configure
                               make
                               make check
                               make install

- Gnuplot
  install from rpm


Install required python modules:
-------------------------------

* Numeric
  install from rpm
  or
  from source
        - download from http://sourceforge.net/projects/numpy
  	- in Numeric-24.0b2:  python setup.py install -- NOTE 1
	
* NetCDF
  install from rpm
  or
  from source
        - download from http://www.unidata.ucar.edu/software/netcdf/
        - in netcdf-3.6.0-p1/src:  ./configure --prefix=/usr/local
                                   make test
                                   make install

* Scientific python
	- download from http://starship.python.net/~hinsen/ScientificPython/
	- in ScientificPython-2.4.9:  python setup.py build
	                              python setup.py install
	
* Biggles
	- download from http://biggles.sourceforge.net/
	- in python2-biggles-1.6.4:  python setup.py build -- NOTE 2
	                             python setup.py install
	
	
* BioPython
        1. mxTextTool
	 - download from http://www.egenix.com/files/python/mxTextTools.html
	 - in egenix-mx-base-2.0.6:  python setup.py install
	2. BioPython
	 - download from http://www.biopython.org/
	 - in biopython-1.40b:  python setup.py build
                                python setup.py test  -- NOTE 3
                                python setup.py install
			   
			   
* PVM
  install pvm and pvm-devel from rpm
  or
  from source
	 - download from http://www.csm.ornl.gov/pvm/pvm_home.html
	 - set the environment variable "PVM_ROOT" to the pvm3 folder
	 - in pvm3: make

  setup environment:
	- example: add to .zshenv:
	export PVM_ROOT=/usr/lib/pvm3  ## depends on the system
	export PVM_RSH=/usr/bin/ssh
	optional: export PVM_ARCH=LINUX 
	         (the variable should match the name of a folder 
		 describing the machine architecture in $PVM_ROOT/lib/)
  
  compile pypvm (C -> Python interface to pvm library):
	  cd biskit/Biskit/PVM
	  make (Note: the Makefile might have to be edited:
	  	VERSION - should matchyour pyhon version
	        LFLAGS - should include to your pvm installation. 
		IFLAGS - the pythom shared libraries might reside in
		         /usr/local/include/python2.4/
			 also check the pvm path
	- optional: compile platform-specific pypvm versions for
	  different python versions or hardware architectures
	  Note: The Biskit.PVM package first tries to import pypvm
		from a subfolder called like the current return value
		of Biskit.tools.platformFolder(). If no such folder
		exists, the default pypvm_core.so is taken from
		biskit/Biskit/PVM/.
	  Example:
	  cd biskit/Biskit/PVM
	  mkdir py2.4_x86_64  (if not already available)
	  cd py2.4_x86_64
	  make
	- see pvm.howto for testing your PVM installation!


* Adapt Biskit settings (binary paths etc):
	- create default settings and adapt them
	  biskit/scripts/Biskit/setup_biskit.py
	  edit ~/.biskit/settings.dat
	

* To use a specific python version (at least 2.4 required):
        - create an alias to the python version that you are 
          going to use (e.g. ln -s /usr/bin/python2.4 python)
        - then edit your enviorment so that your linked 
          python version ends up first in your PATH.
	Example(shell; settings file; export statemant)
	   tcsh ; .~/cshrc   ; setenv PATH ~/name:$PATH
           zsh  ;  ~/.zshenv ; export PATH=~/name:$PATH


* optional: include scripts in search path (zsh, bash)
	BPath=~/biskit/scripts/
	export PATH=$PATH:$BPath/Biskit:$BPath/Dock:$BPath/Mod


Required programs: Biskit
-------------------------

all optional: 
- xplor, gbxplor (compiled for generalized Born)
- HMMER
- whatif
- prosa2003
- foldX
- pymol

Required programs: Biskit/Mod
-----------------------------

- blastall, fastacmd, blastclust (= NCBI blast programs)
- local sequence data bases nr, pdbaa, sprot (formatted with -o)
- local PDB structure data base (folder with *.ent files)
- t_coffee
- modeller

Required programs: Biskit/Dock
------------------------------

- xplor, HEX, WhatIf
- recommended: pymol, prosa, hmmer, foldX


How to use the Mod package:
---------------------------

A test case is included in the biskit CVS. The following should work
without errors (assuming biskit/scripts/Mod is in your search path):

cd ~/biskit/test/Mod/project
search_sequences.py
search_templates.py
clean_templates.py
align.py
model.py

By default, the 5 modelling scripts look for a file "target.fasta" in
the current directory. If that is available, the scripts won't request
any parameters and will create different result folders in the current
folder. Nevertheless, there are a lot of options that can be
specified. Call each script with the option -help to get an overview.

further documentation is pending...


===================== NOTES ===============================

NOTE 1:
======

Lapack is broken in some recent varsions of Numerical Python. This will affect
some Biskit modules (e.g clustering and entropy calculations) as well as some
BioPython modules.
 
To test if your Lapack is broken type this in your python interpreter:
>>> from Numeric import *
>>> from LinearAlgebra import eigenvalues
>>> x = array([[1,2],[3,4]])
>>> eigenvalues(x)

If it hangs you will have to find a version of Numeric without this bug or find find a
way to compile Numeric not resulting in this bug. One way to do this is to
install ALTAL and LAPACK localy and then compile Numeric using this. 
More instructions here: https://cirl.berkeley.edu/view/Grants/BuildLapack

ATLAS: http://math-atlas.sourceforge.net
in ATLAS: make 
          make install arch=<arch> (takes a very long time)

Lapack: http://www.netlib.org/lapack/
in LAPACK:  cp INSTALL/make.inc.LINUX make.inc
            make lapacklib
            make clean
	    
Complementing ATLAS with LAPACK	    
cd ~/ATLAS/lib/<arch>
    cp liblapack.a liblapack_orig.a # make a backup
    mkdir tmp; cd tmp
    ar x ../liblapack.a
    cp ~/LAPACK/lapack_LINUX.a ../liblapack.a
    ar r ../liblapack.a *.o
    cd ..; rm -rf tmp	    

Install:
   mkdir /usr/local/lib/atlas (if not there)
   cp ~/ATLAS/lib/<arch>/*.a  /usr/local/lib/atlas/.
  
   mkdir /usr/local/include/atlas (if not there)
   cp ~/ATLAS/include/cblas.h /usr/local/include/atlas/.
   cp ~/ATLAS/include/clapack.h /usr/local/include/atlas/.


Then in Numeric-24.0b2 edit customize.py:
   use_system_lapack = 1
   lapack_library_dirs = ['/usr/local/lib/atlas']
   lapack_libraries = ['lapack', 'cblas', 'f77blas', 'atlas', 'g2c']
   lapack_extra_link_args = []

Now install Numeric and things should work.
    python setup.py build
    python setup.py install



NOTE 2:
======
Error: 
  /usr/bin/ld: cannot find -lXaw
  (Biggles 1.6.4 on a freshly installed CentOS 4.1 machine)
Cause:
   links in /usr/X11R6/lib missing
Remedy:
   link libXaw.so to libXaw.so.7.0
   (i.e. ln -s /usr/X11R6/lib/libXaw.so.7.0 /usr/X11R6/lib/libXaw.so)
Comment:
   the same link problem occured for libXmu, libXt, libSM, 
   libICE, libXext and libX11


NOTE 3:
======
If the biopython test hags at "test_SVDSuperimposer ..." you are suffering from
the Numeric/eigenvalue problems. See NOTE 1.


How to use the Dock package:
----------------------------

1) Clean up and convert pdb's to xplor psf and pdb
   
   	pdb2xplor 	-i rawRec.pdb     -a -o ./rec_wet/
   	pdb2xplor 	-i rawLig.pdb     -a -o ./lig_wet/ -c 1
   
   	in ./rec_wet/:	gbxplor < generate.inp > xp.log
   	in ./lig_wet/:	-"-
   
   Optional - Reference complex
   	Note: The chain identifiers of the reference complex have to
              match those of the receptor and ligand pdb-files. 
	      The offset of the chain id can be controlled via the -c flag.
	  
   	pdb2xplor 	-i rawComplex.pdb -a -o ./ref_com/   
   	in ./com_wet/:	gbxplor < generate.inp > xp.log


2) Run PCR to generate a number of structures representing the 
   dynamics of the systems.
   
   mkdir rec_pc2_00
   runPcr.py -t ./rec_wet -h host_computer -r ./rec_pc2_00 -f 0
   
   mkdir lig_pc2_00
   runPcr.py -t ./lig_wet -h host_computer -r ./lig_pc2_00 -f 0

   -> runs ensemble MD on host_computer with PCR force constant 0 
      (no restraint) and puts result in ???_pc2_00 

3) Convert frame PDBs into pickled PDBModels (parallized)
 
   in rec_pc2_00:    mkdir model
		     pdbs2model.py -i pcr_00/*pdb -h 30 -o model -a
   in lig_pc2_00:	 -"-

   -> run on 30 hosts, save pickled objects to ./model

4) Create Trajectory object from PDBs or PDBModels

   in rec_pc2_00:    pdb2traj.py -i model/* -r model/????_1_0.model
				 -dr ../rec_wet/????.pdb -f

   in lig_pc2_00:    pdb2traj.py -i model/* -r model/????_1_0.model
				 -dr ../lig_wet/????.pdb -f

   -> pickled Trajectory object with Crystal structure as reference
      ("dry reference", a "wet reference" with identical waters has to 
      be given temporarily with -r). Frames are RMS-fitted to reference

5) prepare folder for docking results
   
	mkdir dock			(better with date dock_DDMM )
	in dock:     mkdir pcr_rec pcr_lig hex01 com

   The Trajectory object  from the dynamics is linked into:   
   	dock/pcr_rec/traj.dat
	dock/pcr_lig/traj.dat

   The following commands assume we are in folder dock.

6) Cluster lig and rec trajectories and select 10 snapshots for docking

   in ./pcr_rec:       selectModels.py -i ./traj.dat -n 10
   in ./pcr_lig:         -"-
   -> 10 pickled PDBModels in pcr_rec and pcr_lig

7) Bundle models into hex-readable pdb file and put them into 
   dictionary of PCRModels.
   
   in ./pcr_rec/: 	PCR2hex   -psf ../../rec_wet/*psf -pdb *model
   in ./pcr_lig/:	  -"-

   in ./com/:		PCR2hex	  -psf ../../com_wet/*psf -pdb complex.pdb
   

8) Create hex macro file
   
   mkdir hex
   in ./hex01/:	hexInput -r ../pcr_rec/rec_hex.pdb
			 -l ../pcr_lig/lig_hex.pdb
			 -c ../com/com_hex.pdb

   edit ./hex01/hex.mac: -docking_sort_mode=0 // for sorting by E
			 -docking_electrostatics=1  // electrostatics on
			 -commenent out: "save_range"  // dont write pdb-files

9) Run hex rigid body docking

   in ./hex01/:	hex -ncpu 2 -nice 0 < hex.mac > hex.log


10) Extract the Hex results from rec-lig_hex.out and save in dictionary

   in ./hex01/:	hex2complex.py 	-rec ../pcr_rec/rec_model.dic 
  				-lig ../pcr_lig/lig_model.dic 
				-hex rec-lig_hex.out 
				-o complexes.dic

11) Calculate protein-protein contact matrices, fnac, scores

   in ./hex01/:	contacter.py -n 10	  // nice level
			     -i *.dic     // complex dic. from hex2complex.py
			     -m contacts  // name of directory to put contact
							matrices in
   
