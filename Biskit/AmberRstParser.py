##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 2 of the
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
## $Revision$
## last $Date$
## last $Author$

"""
Convert Amber restart file to array, PDBModel or Amber crd file.
"""
    
import re
import Numeric as N
import sys
import os.path

from AmberCrdParser import AmberCrdParser, ParseError
import tools as T

class AmberRstParser( AmberCrdParser ):
    """
    Convert Amber restart file to array, PDBModel or Amber crd file.
    """

    def __init__( self, frst, rnAmber=0 ):
        """
        @param frst: input restart file
        @type  frst: str
        @param rnAmber: rename Amber to standard residues (HIE, HID, HIP, CYX)
        @type  rnAmber: 1|0
        """
        self.frst = T.absfile( frst )
        self.crd  = open( self.frst )

        ## number of atoms, number of lines
        self.n = 0
        self.lines_per_frame = 0
        self.xyz = None

        if rnAmber:
            self.renameAmberRes( self.ref )

        ## fields that are not needed by RstParser
        self.ref  = None
        self.box = 0

        ## pre-compile pattern for line2numbers
        xnumber = "-*\d+\.\d+"              # optionally negtive number
        xspace  = ' *'                      # one or more space char
        self.xnumbers = re.compile('('+xspace+xnumber+')')


    def crd2traj( self ):
        """
        @raise ParseError: Not supported for AmberRstParser
        """
        raise ParseError("Not supported for AmberRstParser")


    def getXyz( self ):
        """
        Get coordinate array.
        
        @return: coordinates, N.array( N x 3, 'f')
        @rtype: array
        
        @raise ParseError: if can't interprete second line
        """
        if not self.xyz:

            ## skip first empty line
            self.crd.readline()

            try:
                self.n = self.crd.readline().split()[0]
                self.n = int( self.n )
            except:
                raise ParseError("Can't interprete second line of "+self.frst)

            ## pre-compute lines expected per frame
            self.lines_per_frame = self.n * 3 / 6

            if self.n % 6 != 0:  self.lines_per_frame += 1

            self.xyz = self.nextFrame()

        return self.xyz


    def getModel( self, ref ):
        """
        Get model.
        
        @param ref: reference with same number and order of atoms
        @type  ref: PDBModel
        
        @return: PDBModel
        @rtype: PDBModel
        """
        if not self.xyz:
            self.getXyz()

        result = ref.clone()
        result.setXyz( self.xyz )

        return result


    def getFirstCrdLine( self ):
        """
        Return the first line of Amber crd.
        
        @return: first line of Amber crd formatted coordinate block
        @rtype: str
        """
        if not self.xyz:
            self.getXyz()

        result = ""
        for x in N.ravel( self.xyz )[:10]:
            result += "%8.3f" % x

        return result + "\n"


    def writeCrd( self, fcrd, append=1, lastAtom=None ):
        """
        Write/Append Amber-formatted block of coordinates to a file.
        If a file handle is given, the file will not be closed.
        
        @param fcrd: file to write to
        @type  fcrd: str or file object
        @param append: append to existing file (default: 1)
        @type  append: 0|1
        @param lastAtom: skip all atoms beyond this one (default: None)
        @type  lastAtom: int
        """
        if not self.xyz:
            self.getXyz()

        if type( fcrd ) == file:
            ## take file handle
            f = fcrd
        else:
            ## create new file handle
            new = (mode=='w' or not os.path.exists( T.absfile(fcrd) ))
            mode = 'w'
            if append:
                mode = 'a'
            f = open( T.absfile( fcrd ), mode )
            if new:
                f.write("\n")

        i = 0
        for x in N.ravel( self.xyz ):
            i = i + 1

            f.write( "%8.3f" % x )

            if (i % 10) == 0:
                f.write("\n")

            if lastAtom and i / 3.0 == lastAtom:
                break

        if ((i % 10) != 0):
            f.write("\n")

        if type( fcrd ) != file:
            ## don't close file that was already given
            f.close()


if __name__ == '__main__':

    f='/home/Bis/raik/interfaces/a11/com_pme_ensemble/1b_eq_solvent_voV/sim.rst'

    p = AmberRstParser( f )

    x = p.getXyz()


