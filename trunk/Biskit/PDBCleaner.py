## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2009 Raik Gruenberg & Johan Leckner
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
##
## last $Author$
## last $Date$
## $Revision$

"""
Clean PDB-files so that they can be used for MD.
"""

import Biskit.molUtils as MU
import Biskit.tools as t
from Biskit.PDBModel import PDBModel

import numpy.oldnumeric as N

import copy

class CleanerError( Exception ):
    pass

class PDBCleaner:
    """
    Remove HETAtoms from PDB, replace non-standard AA by closest standard AA.
    Remove non-standard atoms from standard AA residues (and more).
    """

    def __init__( self, fpdb, log=None ):
        """
        @param fpdb: pdb file OR PDBModel
        @type  fpdb: str
        @param log: LogFile object
        @type  log: object
        """
        self.model = PDBModel( fpdb )
        self.log = log


    def logWrite( self, msg, force=1 ):
        if self.log:
            self.log.add( msg )
        else:
            if force:
                print msg

    def remove_multi_occupancies( self ):
        """
        Keep only atoms with alternate A field (well, or no alternate).
        """
        self.logWrite( self.model.pdbCode +
                       ': Removing multiple occupancies of atoms ...')

        i = 0
        to_be_removed = []

        for a in self.model:

            if a['alternate']:
                try:
                    str_id = "%i %s %s %i" % (a['serial_number'], a['name'],
                                              a['residue_name'],
                                              a['residue_number'])

                    if a['alternate'].upper() == 'A':
                        a['alternate'] = ''

                    else:
                        if float( a['occupancy'] ) < 1.0:
                            to_be_removed += [ i ]
                            self.logWrite(
                                'removing %s (%s %s)' % (str_id,a['alternate'],
                                                         a['occupancy']))
                        else:
                            self.logWrite(
                                ('keeping non-A duplicate %s because of 1.0 '+
                                'occupancy') % str_id )

                except:
                    self.logWrite("Error removing duplicate: "+t.lastError() )
            i+=1

        try:
            self.model.remove( to_be_removed )
            self.logWrite('Removed %i atoms' % len( to_be_removed ) )

        except:
            self.logWrite('No atoms with multiple occupancies to remove' )


    def replace_non_standard_AA( self, amber=0, keep=[] ):
        """
        Replace amino acids with none standard names with standard
        amino acids according to L{MU.nonStandardAA}
        
        @param amber: don't rename HID, HIE, HIP, CYX, NME, ACE [0]
        @type  amber: 1||0
        @param keep: names of additional residues to keep
        @type keep:  [ str ]
        """
        standard = MU.atomDic.keys() + keep

        if amber:
            standard.extend( ['HID', 'HIE', 'HIP', 'CYX', 'NME', 'ACE'] )

        replaced = 0

        self.logWrite(self.model.pdbCode +
                      ': Looking for non-standard residue names...')

        resnames = self.model['residue_name']
        for i in self.model.atomRange():

            resname = resnames[i].upper()

            if resname not in standard:
                if resname in MU.nonStandardAA:
                    resnames[i] = MU.nonStandardAA[ resname ]

                    self.logWrite('renamed %s %i to %s' % \
                                  (resname, i, MU.nonStandardAA[ resname ]))
                else:
                    resnames[i] = 'ALA'

                    self.logWrite('Warning: unknown residue name %s %i: ' \
                                  % (resname, i ) )
                    self.logWrite('\t->renamed to ALA.')

                replaced += 1

        self.logWrite('Found %i atoms with non-standard residue names.'% \
                      replaced )


    def __standard_res( self, resname ):
        """
        Check if resname is a standard residue (according to L{MU.atomDic})
        if not return the closest standard residue (according to
        L{MU.nonStandardAA}).
        
        @param resname: 3-letter residue name
        @type  resname: str
        
        @return: name of closest standard residue or resname itself
        @rtype: str
        """
        if resname in MU.atomDic:
            return resname

        if resname in MU.nonStandardAA:
            return MU.nonStandardAA[ resname ]

        return resname


    def remove_non_standard_atoms( self ):
        """
        First missing standard atom triggers removal of standard atoms that
        follow in the standard order. All non-standard atoms are removed too.
        Data about standard atoms are taken from L{MU.atomDic} and symomym
        atom name is defined in L{MU.atomSynonyms}.
        
        @return: number of atoms removed
        @rtype: int
        """
        mask = []

        self.logWrite("Checking content of standard amino-acids...")

        for res in self.model.resList():

            resname  = self.__standard_res( res[0]['residue_name'] ).upper()
            standard = copy.copy( MU.atomDic[ resname ] )

            ## replace known synonyms by standard atom name
            for a in res:
                n = a['name']
                if not n in standard and MU.atomSynonyms.get(n,0) in standard:
                    a['name'] = MU.atomSynonyms[n]
                    self.logWrite('%s: renaming %s to %s in %s %i' %\
                                  ( self.model.pdbCode, n, a['name'],
                                    a['residue_name'], a['residue_number'] ) )

            anames   = [ a['name'] for a in res ]
            keep = 1

            ## kick out all standard atoms that follow a missing one
            rm = []
            for n in standard:
                if not n in anames:
                    keep = 0

                if not keep:
                    rm += [ n ]

            for n in rm:
                standard.remove( n )

            ## keep only atoms that are standard (and not kicked out above)
            for a in res:

                if a['name'] not in standard:
                    mask += [1]
                    self.logWrite('%s: removing atom %s in %s %i '%\
                                  ( self.model.pdbCode, a['name'],
                                    a['residue_name'], a['residue_number'] ) )
                else:
                    mask += [0]

        self.model.remove( mask )

        self.logWrite('Removed ' + str(N.sum(mask)) +
                      ' atoms because they were non-standard' +
                      ' or followed a missing atom.' )

        return N.sum( mask )


    def process( self, keep_hetatoms=0, amber=0, keep_xaa=[] ):
        """
        Remove Hetatoms, waters. Replace non-standard names.
        Remove non-standard atoms.
        
        @param keep_hetatoms: option
        @type  keep_hetatoms: 0||1
        @param amber: don't rename amber residue names (HIE, HID, CYX,..)
        @type  amber: 0||1
        @param keep_xaa: names of non-standard residues to be kept
        @type  keep_xaa: [ str ]
        
        @return: PDBModel (reference to internal)
        @rtype: PDBModel
        
        @raise CleanerError: if something doesn't go as expected ...
        """
        try:
            if not keep_hetatoms:
                self.model.remove( self.model.maskHetatm() )

            self.model.remove( self.model.maskH2O() )

            self.model.remove( self.model.maskH() )

            self.remove_multi_occupancies()

            self.replace_non_standard_AA( amber=amber, keep=keep_xaa )

            self.remove_non_standard_atoms()


        except KeyboardInterrupt, why:
            raise KeyboardInterrupt( why )
        except Exception, why:
            self.logWrite('Error: '+t.lastErrorTrace())
            raise CleanerError( 'Error cleaning model: %r' % why )

        return self.model



#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test class """

    def prepare(self):
        from Biskit.LogFile import LogFile
        import tempfile

        self.f_out = tempfile.mktemp( '_test_PDBCleaner' )
        
        self.l = LogFile( self.f_out, mode='w')


    def test_PDBCleaner( self ):
        """PDBCleaner test"""
        
        ## Loading PDB...
        self.c = PDBCleaner( t.testRoot() + '/rec/1A2P_rec_original.pdb',
                             log=self.l )
        
        self.m = self.c.process()

        if self.local:
            print 'PDBCleaner log file written to: %s'%self.f_out

        self.assertAlmostEqual( N.sum( self.m.masses()), 34029.0115499993, 7 )

    def cleanUp(self):
        t.tryRemove( self.f_out )    
    
        
if __name__ == '__main__':

    BT.localTest()
