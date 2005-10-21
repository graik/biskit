(C) Oleg Tsodikov, 2001

Surface Racer 3.0 is distributed free of charge only for academic use. 

System requirements: Windows 95/98/NT/2000


INSTALLATON


WINDOWS
After the download, use the EXtract All  command (in the WinZip utility, which should start
automatically after clicking on sracer.zip) to extract fastsurf.exe, surfrace.exe, radii.txt, 
and readme.txt from the archive sracer.zip 
into the same directory. The programs (fastsurf.exe and surfrace.exe) are now ready 
to be used. 

UNIX-like systems (LINUX and OSX)
Allow the .exe files to be executed (by the chmod command) and they are ready to run.

INPUT AND OUTPUT

Surface Racer  uses files in the Protein Data Bank format as the input. There is no
limit on the number of atoms in the structure, but side chains (no hydrogen atoms should be present!) 
need to be included. 
The output file has the same name with the ".txt" suffix and contains only ATOM records 
in the original format, with van der Waals radius, accessible surface area, molecular 
surface area, average curvature of molecular surface for each atom listed at the end of 
each line in this particular order, the curvature being the last entry. The same information
for the cavities in protein interior is recorded at the end of the file (when using surfrace). Note, that sometimes
the same atom can be accessible from both outside and a cavity if this cavity is 
sufficiently close to the surface.  

In addition, the file result.txt is created. This file contains the probe radius, total accessible
and molecular surface areas and their breakdown into polar, nonpolar, charged etc parts.
These results are added to the file after each new calculation.

File radii.txt contains two van der Waals radus sets, commonly used in macromolecular 
surface area calculations. The radius values can be changed by the user provided that the
file format is fully preserved.

You can perform calculation of accessible and molecular surface area with either fastsurf
and surfrace. Surfrace has the additional feature of cavity analysis, but because of that
it runs a little more slowly.

In case of a file containing two or more disjointed structures, only the topmost structure (the molecule containing the atom
with the largest z coordinate) will be calculated properly therefore it is best to place these structures
in individual pdb files and run the program sepapately for each file.

CALCULATION

The program runs as a console application. The user chooses between the first and the
second radus sets from radii.txt (by pressing 1 or 2), inputs the name of the PDB file (the
input file should be in the same directory as the program and should include the suffix), 
probe radius in angstroms and the calculation mode. The calculation is performed only
for atoms with the ATOM record. If the file contains disconnected structures (where "disconnected" 
means that they cannot be bridged by the probe sphere), only one of these structures will 
be included in the calculation. The safest way in this case is to leave only the structure of 
interest in the file before running Surface Racer.

REFERENCE

The algorithm employed in this program has been described in detail in the following 
paper:

Tsodikov, O. V., Record, M. T. Jr. and Sergeev, Y. V. (2002). A novel computer program 
for fast exact calculation of accessible and molecular surface areas and average surface 
curvature. J. Comput. Chem., 23, 600-609. 


