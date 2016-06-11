## Convert PDB files to a pickled PDBModel object.
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2016 Raik Gruenberg & Johan Leckner
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
##
## required by pdbs2struct.py

"""
Convert a PDB file to a pickled PDBModel object.
"""
    
from PDBModel import PDBModel
import tools as T

## PVM imports
from Biskit.PVM import JobSlave


class StructureSlave(JobSlave):
    """
    Convert a PDB file to a pickled PDBModel object.
    """

    def initialize(self, params):
        """
        Copy the parameters that Master is passing in as dict into
        fields of this class.
        
        @param params: defined in Master
        @type  params: dict
        """
        self.params = params


    def renameAmberRes( self, model ):
        """
        Rename special residues (from Amber) back into standard
        names (i.e CYX S{->} CYS).
        
        @param model: model
        @type  model: PDBModel
        """
        for a in model.getAtoms():
            if a['residue_name'] == 'CYX':
                a['residue_name'] = 'CYS'
            if a['residue_name'] in ['HIE','HID','HIP']:
                a['residue_name'] = 'HIS'


    def renameToAmberAtoms( self, model ):
        """
        ptraj puts last letter/number of 4-letter atom names first. Undo.
        
        @param model: model
        @type  model: PDBModel
        
        @note: could be avoided if ptraj would be told:
               trajout |file| pdb nowrap
        """
        numbers = map( str, range(10) )

        for a in model.getAtoms():
            if len(a['name'])==4 and a['name'][0] in numbers:
                a['name'] = a['name'][1:] + a['name'][0]


    def go(self, dict):
        """
        Calculate rmsd values for trajectory and plot them.

        @param dict: dictionary with path to pdb files as keys
        @type  dict: dict

        @return: dictionary mapping input:output names
        @rtype: dict        
        """
        d = {}

        print "working on :",
        for pdbIn, out in dict.items():

            print T.stripFilename( pdbIn )

            ## set file name for pickled Structure object
            dir = os.path.dirname( pdbIn )
            if self.params['out'] <> None:
                dir = self.params['out']

            out = dir + '/' + T.stripFilename( pdbIn) + '.model'

            try:
                m = PDBModel( pdbIn, skipRes=self.params['skipRes'] )

                if self.params['amber']:
                    self.renameAmberRes( m )
                    self.renameToAmberAtoms( m )

                if self.params['sort']:
                    m = m.sort()

                m.saveAs( out )

            except:
                if self.params.has_key('report') and self.params['report']:
                    f = open( out[:-6] + '.error', 'w' )
                    f.write( T.lastError() )
                    f.close()
                else:
                    T.errWriteln("Error" + T.lastError()  )


                out = ''

            d[pdbIn] = out

        print "Done."
        return d

##############
## empty test
##############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Mock test, the Slave is tested in L{Biskit.StructureMaster}"""
    pass

if __name__ == '__main__':

    import os, sys

    if len(sys.argv) == 2:

        niceness = int(sys.argv[1])
        os.nice(niceness)

    slave = StructureSlave()
    slave.start()
