## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py
## DAG - substituted Numeric

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2012 Raik Gruenberg & Johan Leckner
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
## $Revision$
## last $Date$
## last $Author$
"""
Convert single amber crd into Trajectory object
"""

import re
import numpy as N
import sys

import tools as T
from Trajectory import Trajectory
from PDBModel import PDBModel
from LogFile import StdLog

def _use():

    print """
Convert single amber crd into Trajectory object

amber2traj.py -i sim.crd -o traj_0.dat -r ref.pdb [-b -wat -hyd -rnres
              -code PDBC ]

    -b     traj has box info (3 additional coordinates per frame)
    -wat   delete WAT, Cl-, Na+ residues (after parsing)
    -hyd   delete all hydrogens (after parsing)
    -rnres rename amber residues HIE/HID/HIP, CYX to HIS and CYS
    -code  PDB code of molecule (otherwise first 4 letters of ref file name)

    ref.pdb must have identical atom content as sim.crd
    """
    sys.exit( 0 )

class ParseError( Exception ):
    pass


class AmberCrdParser:
    """
    Convert an Amber-generated crd file into a Trajectory object.
    """

    def __init__( self, fcrd, fref, box=0, rnAmber=0, pdbCode=None,
                  log=StdLog(), verbose=0 ):
        """
        @param fcrd: path to input coordinate file
        @type  fcrd: str
        @param fref: PDB or pickled PDBModel with same atom content and order
        @type  fref: str
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
        self.fcrd = T.absfile( fcrd )
        self.crd  = T.gzopen( self.fcrd )

        self.ref  = PDBModel( T.absfile(fref), pdbCode=pdbCode )
        self.box  = box

        self.n = self.ref.lenAtoms()

        self.log = log
        self.verbose = verbose

        if rnAmber:
            self.ref.renameAmberRes()

        ## pre-compile pattern for line2numbers
        xnumber = "-*\d+\.\d+"              # optionally negtive number
        xspace  = ' *'                      # one or more space char
        self.xnumbers = re.compile('('+xspace+xnumber+')')

        ## pre-compute lines expected per frame
        self.lines_per_frame = self.n * 3 / 10

        if self.n % 10 != 0:  self.lines_per_frame += 1
        if self.box:          self.lines_per_frame += 1

        ## mark chains (the TER position might get lost when deleting atoms)
        ## should not be necessary any longer
##        if not self.ref.getAtoms()[0].get('chain_id',''):
##            self.ref.addChainId()

    def line2numbers( self, l ):
        """
        convert a line from crd/vel file to list of float numbers
        
        @param l: line
        @type  l: str
        
        @return: list of floats
        @rtype: [float]
        """
        match = self.xnumbers.findall( l )

        return [ round( float(strCrd),3) for strCrd in match ] 


    def nextLine( self ):
        """
        extract next 10 coordinates from crd file

        @return: coordinates
        @rtype: [float]    
        """
        l = self.crd.readline()
        if l == '':
            raise EOFError('EOF')

        return self.line2numbers( l )


    def nextFrame( self ):
        """
        Collect next complete coordinate frame

        @return: coordinate frame
        @rtype: array
        """

        i = 0
        xyz = []
        while i != self.lines_per_frame:

            if self.box and i+1 == self.lines_per_frame:

                ## skip box info
                if len( self.nextLine() ) != 3:
                    raise ParseError( "BoxInfo must consist of 3 numbers." )

            else:
                xyz += self.nextLine()

            i += 1

        return N.reshape( xyz, ( len(xyz) / 3, 3 ) ).astype(N.float32)


    def crd2traj( self ):
        """
        Convert coordinates into a Trajectory object.

        @return: trajectory object
        @rtype: Trajectory
        """
        ## skip first empty line
        self.crd.readline()

        xyz = []
        i = 0

        if self.verbose: self.log.write( "Reading frames .." )

        try:
            while 1==1:

                xyz += [ self.nextFrame() ]
                i += 1

                if i % 100 == 0 and self.verbose:
                    self.log.write( '#' )

        except EOFError:
            if self.verbose: self.log.add("Read %i frames." % i)

        t = Trajectory( refpdb=self.ref )

        t.frames = N.array( xyz ).astype(N.float32)

        t.setRef( self.ref )
        t.ref.disconnect()

        return t

import Biskit.test as BT
import tempfile

class Test( BT.BiskitTest ):
    """Test AmberCrdParser"""

    def prepare(self):
        root = T.testRoot() + '/amber/'
        self.finp = root + 'sim.crd.gz'
        self.fref = root + '1HPT_0.pdb'
        self.fout = tempfile.mktemp('.dat', 'traj_')

    def cleanUp(self):
        T.tryRemove( self.fout )
        

    def test_AmberCrdParser(self):
        """AmberCrdParser test"""
        
        self.p = AmberCrdParser( self.finp, self.fref, box=True, rnAmber=True,
                                 log=self.log, verbose=self.local )
        self.t = self.p.crd2traj()

        self.t.removeAtoms(lambda a: a['residue_name'] in ['WAT','Na+','Cl-'] )
        self.t.removeAtoms(lambda a: a['element'] == 'H' )

        if self.local:
            print "Dumping result to ", self.fout

        T.dump( self.t, T.absfile(self.fout) )

        if self.local:
            print "Dumped Trajectory with %i frames and %i atoms." % \
                  (len(self.t), self.t.lenAtoms() )

        self.assertEqual( len(self.t), 10 )
        self.assertEqual( self.t.lenAtoms(), 440 )

if __name__ == '__main__':

    BT.localTest()
