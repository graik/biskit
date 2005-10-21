(C) Oleg Tsodikov, 2001

Surface Racer 1.1 is distributed free of charge only for academic use. 

System requirements: Unix/Linux


INSTALLATON

Make sure that readme, radii, surfrace, and fastsurf files are downloaded into the same
directory. Type: chmod +x surfrace (or some other analogous command) to allow surfrace 
to be executed. Repeat the same thing for fastsurf. Run programs by typing ./surfrace
or ./fastsurf.

INPUT AND OUTPUT

Surface Racer 1.1 uses files in the Protein Data Bank format as the input. There is no
limit on the number of atoms in the structure, but side chains need to be included. 
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

CALCULATION

The user chooses between the first and the
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


