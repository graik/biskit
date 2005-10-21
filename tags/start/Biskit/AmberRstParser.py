##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## $Revision$
## last $Date$
## last $Author$

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
        frst - str, input restart file
        fref - str, PDB or pickled PDBModel with same atom content and order
        rnAmber - 1|0, rename Amber to standard residues (HIE, HID, HIP, CYX)
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
        raise ParseError("Not supported for AmberRstParser")


    def getXyz( self ):
        """
        -> N.array( N x 3, 'f')
        !! ParseError
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
        ref - PDBModel, reference with same number and order of atoms 
        -> PDBModel
        """
        if not self.xyz:
            self.getXyz()

        result = ref.clone()
        result.setXyz( self.xyz )

        return result


    def getFirstCrdLine( self ):
        """
        -> str, first line of Amber crd formatted coordinate block
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
        fcrd   - str or file object, file to write to
        append - 0|1, append to existing file [1]
        lastAtom - int, skip all atoms beyond this one [None]
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

    
