## Convert PDB files to a pickled PDBModel object.
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## required by pdbs2struct.py
##
## $Revision$
## last $Date$
## last $Author$

from PDBModel import PDBModel
import tools as T

## PVM imports
from PVM.dispatcher import JobSlave


class StructureSlave(JobSlave):
    """
    Convert a PDB file to a pickled PDBModel object.
    """

    def initialize(self, params):

        self.params = params


    def renameAmberRes( self, model ):

        for a in model.getAtoms():
            if a['residue_name'] == 'CYX':
                a['residue_name'] = 'CYS'
            if a['residue_name'] in ['HIE','HID','HIP']:
                a['residue_name'] = 'HIS'


    def renameToAmberAtoms( self, model ):
        """ptraj puts last letter/number of 4-letter atom names first. Undo.
        could be avoided if ptraj would be told: trajout |file| pdb nowrap
        """
        numbers = map( str, range(10) )

        for a in model.getAtoms():
            if len(a['name'])==4 and a['name'][0] in numbers:
                a['name'] = a['name'][1:] + a['name'][0]


    def go(self, dict):

        d = {}

        print "working on :",
        for pdbIn, out in dict.items():

            print T.stripFilename( pdbIn )
            
            ## set file name for pickled Structure object
            dir = os.path.dirname( pdbIn)
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
    
if __name__ == '__main__':

    import os, sys

    if len(sys.argv) == 2:
        
        niceness = int(sys.argv[1])
        os.nice(niceness)

    slave = StructureSlave()
    slave.start()
