## numpy-oldnumeric calls replaced by custom script; 09/06/2016
## Automatically adapted for numpy-oldnumeric Mar 26, 2007 by alter_code1.py

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2018 Raik Gruenberg & Johan Leckner
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
"""
Parse Amber restart files.
"""

import re
import os.path

import biskit.core.oldnumeric as N0
from biskit import PDBModel
import biskit.tools as T

## allow relative imports when calling module by itself for testing (pep-0366)
if __name__ == "__main__" and __package__ is None:
    import biskit.md; __package__ = "biskit.md"

class ParseError(Exception):
    pass

class AmberRstParser:
    """Convert an Amber restart file to array, PDBModel or a Amber crd file.

    Note: AmberRstParser is currently ignoring both the velocity and
    boxinfo record (although this could be easily changed).
    """

    def __init__( self, frst ):
        """
        :param frst: input restart file
        :type  frst: str
        """
        self.frst = T.absfile( frst )
        self.crd  = open( self.frst )

        self.n = 0       #: number of atoms
        self.lines_per_frame = 0
        self.xyz = None  #: will hold coordinate array
        self.box = None  #: will hold box array if any

        ## pre-compile pattern for line2numbers
        xnumber = r"-*\d+\.\d+"              # optionally negtive number
        xspace  = r' *'                      # one or more space char
        self.xnumbers = re.compile('('+xspace+xnumber+')')


    def __del__(self):
        try:
            self.crd.close()
        except:
            pass


    def __nextLine( self ):
        """Extract next line of coordinates from crd file

        :return: coordinates
        :rtype: [float]    
        """
        l = self.crd.readline()
        if l == '':
            raise EOFError('EOF')

        match = self.xnumbers.findall( l )
        return [ round( float(strCrd),7) for strCrd in match ] 


    def __frame( self ):
        """Collect next complete coordinate frame

        :return: coordinate frame
        :rtype: array
        """
        self.xyz = [ self.__nextLine() for i in range(self.lines_per_frame) ]

        return N0.reshape(self.xyz, ( self.n, 3 ) ).astype(N0.Float32)


    def getXyz( self ):
        """Get coordinate array.

        :return: coordinates, N0.array( N x 3, 'f')
        :rtype: array

        :raise ParseError: if can't interprete second line
        """
        if not self.xyz:

            ## skip first empty line
            self.crd.readline()

            try:
                self.n, self.time = self.crd.readline().split()
                self.n = int( self.n )
                self.time = float( self.time )
            except:
                raise ParseError("Can't interprete second line of "+self.frst)

            ## pre-compute lines expected per frame
            self.lines_per_frame = self.n // 2
            if self.n % 2 != 0:
                self.lines_per_frame += 1

            self.xyz = self.__frame()

        return self.xyz


    def getModel( self, ref, rnAmber=0 ):
        """
        Get model.

        :param ref: reference with same number and order of atoms
        :type  ref: PDBModel
        :param rnAmber: rename Amber to standard residues (HIE, HID, HIP, CYX)
        :type  rnAmber: 1|0

        :return: PDBModel
        :rtype: PDBModel
        """
        if not self.xyz:
            self.getXyz()

        result = ref.clone()
        result.setXyz( self.xyz )
        if rnAmber:
            result.renameAmberRes()

        return result


    def getFirstCrdLine( self ):
        """
        Return the first line of Amber crd.

        :return: first line of Amber crd formatted coordinate block
        :rtype: str
        """
        if not self.xyz:
            self.getXyz()

        result = ""
        for x in N0.ravel( self.xyz )[:10]:
            result += "%8.3f" % x

        return result + "\n"


    def writeCrd( self, fcrd, append=1, lastAtom=None ):
        """
        Write/Append Amber-formatted block of coordinates to a file.
        If a file handle is given, the file will not be closed.

        :param fcrd: file to write to
        :type  fcrd: str or file object
        :param append: append to existing file (default: 1)
        :type  append: 0|1
        :param lastAtom: skip all atoms beyond this one (default: None)
        :type  lastAtom: int
        """
        if not self.xyz:
            self.getXyz()

        if type( fcrd ) == file:
            ## take file handle
            f = fcrd
        else:
            ## create new file handle
            mode = 'w'
            if append:
                mode = 'a'
            f = open( T.absfile( fcrd ), mode )

            newf = (mode=='w' or not os.path.exists( T.absfile(fcrd) ))
            if newf:
                f.write("\n")

        i = 0
        for x in N0.ravel( self.xyz ):
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

######################
### Module testing ###
import biskit.test as BT

class Test(BT.BiskitTest):
    """Test AmberRstParser"""

    def prepare(self):
        self.f    = T.testRoot()+'/amber/sim.rst'
        self.fref = T.testRoot()+'/amber/1HPT_0.pdb'

        self.p = AmberRstParser( self.f )

    def test_getXyz( self ):
        """AmberRstParser.getXyz test"""
        self.xyz = self.p.getXyz()
        self.assertEqual( N0.shape(self.xyz), (11200,3) )

class TestLong(Test):
    """long AmberRstParser test"""
    TAGS = [BT.LONG]

    def test_getModel(self):
        """AmberRstParser.getModel test"""
        self.ref = PDBModel( self.fref )
        self.model = self.p.getModel( self.ref )
        self.assertEqual( len(self.model), 11200 )


if __name__ == '__main__':

    ## run Test and push self.* fields into global namespace
    BT.localTest( )


