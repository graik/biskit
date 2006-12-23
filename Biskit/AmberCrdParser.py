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
Convert single amber crd into Trajectory object
"""

import re
import Numeric as N
import sys

import tools as T
from Trajectory import Trajectory
from PDBModel import PDBModel

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

    def __init__( self, fcrd, fref, box=0, rnAmber=0, pdbCode=None ):
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
        """
        self.fcrd = T.absfile( fcrd )
        self.crd  = open( self.fcrd )

        self.ref  = PDBModel( T.absfile(fref), pdbCode=pdbCode )
        self.box  = box

        self.n = self.ref.lenAtoms()

        if rnAmber:
            self.renameAmberRes( self.ref )

        ## pre-compile pattern for line2numbers
        xnumber = "-*\d+\.\d+"              # optionally negtive number
        xspace  = ' *'                      # one or more space char
        self.xnumbers = re.compile('('+xspace+xnumber+')')

        ## pre-compute lines expected per frame
        self.lines_per_frame = self.n * 3 / 10

        if self.n % 10 != 0:  self.lines_per_frame += 1
        if self.box:          self.lines_per_frame += 1

        ## mark chains (the TER position might get lost when deleting atoms)
        if not self.ref.getAtoms()[0].get('chain_id',''):
            self.ref.addChainId()


    def renameAmberRes( self, model ):
        """
        Rename special residue names from Amber back into standard names
        (i.e CYX S{->} CYS )
        
        @param model: model, will be modified in-place
        @type  model: PDBModel
        """
        for a in model.getAtoms():
            if a['residue_name'] == 'CYX':
                a['residue_name'] = 'CYS'
            if a['residue_name'] in ['HIE','HID','HIP']:
                a['residue_name'] = 'HIS'


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

        return N.reshape( xyz, ( len(xyz) / 3, 3 ) ).astype(N.Float32)


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

        T.flushPrint( "Reading frames .." )

        try:
            while 1==1:

                xyz += [ self.nextFrame() ]
                i += 1

                if i % 100 == 0:
                    T.flushPrint( '#' )

        except EOFError:
            print "Read %i frames." % i

        t = Trajectory( refpdb=self.ref )

        t.frames = N.array( xyz ).astype(N.Float32)

        t.setRef( self.ref )
        t.ref.disconnect()

        return t


if __name__ == '__main__':

    if len( sys.argv ) < 2:
        _use()

    o = T.cmdDict( {'o':'traj_0.dat', 'i':'sim.crd'} )
    fcrd = o['i']
    fpdb = o['r']
    fout = o['o']
    box  = o.has_key( 'b' )
    wat  = o.has_key('wat')
    hyd  = o.has_key('hyd')
    rnres  = o.has_key('rnres')
    code = o.get('code', None)

    p = AmberCrdParser( fcrd, fpdb, box, rnres, pdbCode=code )
    t = p.crd2traj()

    if wat:
        t.removeAtoms( lambda a: a['residue_name'] in ['WAT', 'Na+', 'Cl-'] )

    if hyd:
        t.removeAtoms( lambda a: a['element'] == 'H' )

    print "Dumping result to ", fout
    T.Dump( t, T.absfile(fout) )

    print "Done"
