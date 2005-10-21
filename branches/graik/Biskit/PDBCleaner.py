##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
##
## last $Author$
## last $Date$
## $Revision$

import Biskit.molUtils as MU
import Biskit.tools as t
from Biskit.PDBModel import PDBModel

from Numeric import sum
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
        fpdb  - str, pdb file OR PDBModel
        log   - LogFile object
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

        for a in self.model.getAtoms():

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


    def replace_non_standard_AA( self, amber=0 ):
        """
        amber - 1|0, don't rename HID, HIE, HIP, CYX, NME, ACE [0]
        """
        standard = MU.atomDic.keys()

        if amber:
            standard.extend( ['HID', 'HIE', 'HIP', 'CYX', 'NME', 'ACE'] )

        replaced = 0

        self.logWrite(self.model.pdbCode +
                      ': Looking for non-standard residue names...')
        
        for a in self.model.getAtoms():

            resname = a['residue_name'].upper()

            if resname not in standard:
                if resname in MU.nonStandardAA:
                    a['residue_name'] = MU.nonStandardAA[ resname ]

                    self.logWrite('renamed %s %i to %s' % \
                                  (resname, a['residue_number'],
                                   MU.nonStandardAA[ resname ]))
                else:
                    a['residue_name'] = 'ALA'

                    self.logWrite('Warning: unknown residue name %s %i: ' \
                                  % (resname, a['residue_number'] ) )
                    self.logWrite('\t->renamed to ALA.')
                    
                replaced += 1

        self.logWrite('Found %i atoms with non-standard residue names.'% \
                      replaced )


    def __standard_res( self, resname ):
        """
        resname - str, 3-letter residue name
        -> str, name of closest standard residue or resname itself
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
        -> int, number of atoms removed
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

        self.logWrite('Removed ' + str(sum(mask)) +
                      ' atoms because they were non-standard' +
                      ' or followed a missing atom.' )

        return sum( mask )
                    

    def process( self, keep_hetatoms=0, amber=0 ):
        """
        Remove Hetatoms, waters. Replace non-standard names.
        Remove non-standard atoms.
        keep_hetatoms - 0|1
        amber         - 0|1, don't rename amber residue names (HIE,HID,CYX,..)
        -> PDBModel (reference to internal)
        !! CleanerError
        """

        try:
            if not keep_hetatoms:
                self.model.remove( self.model.maskHetatm() )

            self.model.remove( self.model.maskH2O() )

            self.model.remove( self.model.maskH() )

            self.remove_multi_occupancies()

            self.replace_non_standard_AA( amber=amber )

            self.remove_non_standard_atoms()


        except KeyboardInterrupt, why:
            raise KeyboardInterrupt( why )
        except Exception, why:
            raise CleanerError( t.lastError() )

        return self.model


if __name__ == '__main__':
   
    c = PDBCleaner( t.testRoot() + '/com/1BGS_original.pdb'  )

    m = c.process()
    
