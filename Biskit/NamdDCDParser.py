## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2009 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.
##
## You find a copy of the GNU General Public License in the file
## license.txt along with this program; if not, write to the Free
## Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
##
##
## $Revision: 742 $
## last $Date: 2009-05-09 14:17:28 +0200 (Sat, 09 May 2009) $
## last $Author: graik $
"""
Convert single amber crd into Trajectory object
"""

import os
import sys
import struct

import numpy as npy
import tools as T
from PDBModel import PDBModel
from LogFile import StdLog

class NamdDCDParser:
    
    def __init__( self, fdcd, fref, box=0, pdbCode=None,
                  log=StdLog(), verbose=0):
        """
        @param fdcd: path to input dcd file
        @type  fdcd: str
        @param fref: PDB or pickled PDBModel or directly an open PDBModel instancewith same atom content and order
        @type  fref: str or PDBModel
        @param box: expect line with box info at the end of each frame
                    (default: 0)
        @type  box: 1|0
        @param rnAmber: rename amber style residues into standard (default: 0)
        @type  rnAmber: 1|0
        @param pdbCode: pdb code to be put into the model (default: None)
        @type  pdbCode: str
        @param log: LogFile instance [Biskit.StdLog]
        @type  log: Biskit.LogFile
        @param verbose: print progress to log [0]
        @type  verbose: int
        """
        self.fdcd = T.absfile( fdcd )
        self.dcd = open(self.fdcd, "r", 0)
        
        if isinstance(fref, str) :
            self.ref=PDBModel(T.absfile(fref), pdbCode=pdbCode)
        elif fref :
            self.ref = fref
            
        self.box  = box
        self.n = self.ref.lenAtoms()
        self.log = log
        self.verbose = verbose
        
        self.readHeader()
        self.set_pointerInfo()

    def readHeader(self, verbose=False):
        """
        Read NAMD DCD coordinates file
        Only for 32bit DCD with opposite endianness!
        AND only if no atom is fixed!
        """
        f = self.dcd
        unpack = struct.unpack
        
        # Check we are in the beggining of the file
        # or move there to read the header
        if f.tell() != 0:
            f.seek(0)
        
        # First read header information and check correct file format
        header = struct.unpack(">I 4s 9I f 11I", f.read(92))
        if header[0] == 84 and header[1] == 'CORD' and header[-1] == 84 and header[-2] != 0:
            #print "recognized 32 bit DCD file of opposite endianness"
            #Store the number of sets of coordinates (NSET). Frames
            self.nset = header[2]
            #Store ISTART, the starting timestep
            self.istart = header[3]
            #Store NSAVC, the number of steps between dcd saves
            self.nsavc = header[4]
            # Store NTOT, number of total simulation steps
            self.ntot = header[5]
            #Store DELTA, the time step of simulation
            self.delta = header[11]
            # Have box information?
            self.has_extrablock = bool(header[12])
            # Have a 4th dimension?
            self.had_4dims = bool(header[13])
        else:
            f.close()
            sys.exit("Bad DCD Format")
            
        # Read title information
        if (unpack('>I',f.read(4))[0] - 4) % 80 == 0:
            # Number of title lines
            self.ntitle = int(unpack('>I', f.read(4))[0])
            self.title = [unpack(">80s",f.read(80)) for i in range(self.ntitle)]
            f.read(4) # Skip closing block number
        else:
            f.close()
            sys.exit("ERROR in title. Bad DCD format")
        
        # Read number of atoms
        atomBlock = unpack('>3I', f.read(12))
        if atomBlock[-1] == 4:
            self.natoms = atomBlock[1]
        else:
            f.close()
            sys.exit("Bad DCD format")
        
        if verbose:
            print self.title
            print "Number of atoms:", self.natoms
            print "Number of frames:",self.nset
            print "Starting timestep", self.istart
            print "Final timestep:", self.ntot
            print "Steps between frames:", self.nsavc
            print "Time step of simulation:", self.delta
    
    def set_pointerInfo(self):
        """
        Store sizes for browsing the file later
        """
        # Header size is: 116 + 80* self.ntitle
        self.h_size = 116 + (80 * self.ntitle)
        
        # Frame size
        # 4 bytes because it's floats per 3 axis per total num of atoms
        # Add the enclosing integers (two for each axis) = 6 * 4
        f_size = (3 * 4 * self.natoms) + 24
        if self.has_extrablock:
            f_size += 56
        
        self.f_size = f_size
           
    def read_charmm_extrablock(self):
        
        f = self.dcd
        unpack = struct.unpack
        
        # This block contains the box information
        if unpack('>I', f.read(4))[0] == 48:
            self.unitcell = npy.fromstring( f.read(48), dtype=">d")
            f.read(4)
        else:
            f.close()
            sys.exit("ERROR in read_charmm_extrablock(). Bad DCD Format")
            
    def read_dcdstep(self):
      
        f = self.dcd
        size = struct.calcsize
        
        # If there is box information
        if self.has_extrablock:
            self.read_charmm_extrablock()
        
        # Read coordinates
        # Each coordinates block is enclosed by one integer
        # that we will skip all the times
        xyz = npy.zeros([self.natoms, 3], dtype=">f4")
        f.read(4)
        xyz[:,0] = npy.fromstring(f.read(size('f')*self.natoms), dtype=">f4")
        f.read(8)
        xyz[:,1] = npy.fromstring(f.read(size('f')*self.natoms), dtype=">f4")
        f.read(8)
        xyz[:,2] = npy.fromstring(f.read(size('f')*self.natoms), dtype=">f4")
        f.read(4)
        
        return xyz
        
    def read_all(self):
        """
        Read all snapshots
        """
        # Go to the beggining of the frames
        f = self.dcd
        f.seek(self.h_size)
        
        # Read Frames
        all_snap = npy.zeros([self.nset, self.natoms, 3], dtype=">f4")
        for i in range(self.nset):
            all_snap[i,:] = self.read_dcdstep()
        
        return all_snap
    
    def close(self):
        self.dcd.close()
        
    def getFrame(self, i):
        """
        Read specific frame.
                """
        f = self.dcd

        # Calculate pointer position
        # for the desired frame and move there
        pointer = self.h_size + (self.f_size * i)
        f.seek(pointer)
        
        return self.read_dcdstep()
    
    def __getitem__(self, i):
        return self.getFrame(i)
 

import Biskit.test as BT
import tempfile

class Test( BT.BiskitTest ):
    #"""Test NamdDCDParser"""
    
    def __init__(self):
        self.pdb = PDBModel("free/HEWL_ETOH_FREE_ETA_0.pdb")
        self.dcd = NamdDCDParser("/home/daniel/WORK/SOLVATION/MDM2/ETA/md/md_5.dcd",self.pdb,box=1)
        self.dcd.readHeader(verbose=1)
        print "Frame 25:",self.dcd.getFrame(24)
        print "Box fr.25:",self.dcd.unitcell
        print "Frame 1:", self.dcd[0]
        print "Box fr.1:",self.dcd.unitcell
        print "Next Frame (2):", self.dcd.read_dcdstep()
        print "all frames:", self.dcd.read_all()
        self.dcd.close()


if __name__ == '__main__':
    
    Test()
    
