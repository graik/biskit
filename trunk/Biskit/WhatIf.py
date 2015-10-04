## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py

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
Run WhatIf, calculate and collect accessibility data
"""

import numpy as N
import string
import tempfile

from Biskit import Executor, TemplateError
## import Biskit.settings as S
import Biskit.tools as T


class WhatIf_Error( Exception ):
    pass

class WhatIf( Executor ):
    """
    WhatIf
    ======
    Run WhatIf accessibility calculation
    
    Result
    ------
    A WhatIf log file containing the accessibilities (in square Angstrom)
    per residue. The ASA is specified for the entire residue as well for
    the backbone and side chain.

    Reference
    ---------
    U{ http://swift.cmbi.kun.nl/whatif/}    

    @note: WhatIf is not free software, therefore the calculations performed
           by this class can now also be performed by the freeware program
           SurfaceRacer, see L{ Biskit.SurfaceRacer }. 
    """

    whatif_script ="""SOUP  \nINISOU

GETMOL  \n%(f_pdb)s  \nTEST

ACCESS  \nINIACC
SETACC  \nALL  \n0   \n     

DOLOG  \n%(f_relativeASA)s  \n0
VACACC  \nN  \nALL
NOLOG

DOLOG  \n%(f_residueASA)s  \n0
SHOACC  \nALL  \n0
NOLOG""" 


    def __init__( self, model, **kw ):
        """
        @param model: PDBModel
        @type  model: 

        @param kw: additional key=value parameters for Executor:
        @type  kw: key=value pairs
        ::
          debug    - 0|1, keep all temporary files (default: 0)
          verbose  - 0|1, print progress messages to log (log != STDOUT)
          node     - str, host for calculation (None->local) NOT TESTED
                          (default: None)
          nice     - int, nice level (default: 0)
          log      - Biskit.LogFile, program log (None->STOUT) (default: None)
        """
        Executor.__init__( self, 'whatif', template=self.whatif_script,
                           f_out='/dev/null', **kw )

        self.f_pdb = tempfile.mktemp('_whatif.pdb')
        self.f_relativeASA = tempfile.mktemp('_whatif_relative.log')
        self.f_residueASA = tempfile.mktemp('_whatif_residue.log')

        self.model = model.clone()


    def prepare( self ):
        """
        Overrides Executor method.
        """
        self.__prepareModel( self.model, self.f_pdb ) 


    def __prepareModel( self, m, f_pdb_out ):
        """
        Prepare a model that WhatIf likes.
         - consecutive numbering of residues
         - remove hydrogens (will be deleted anyway)
        
        @param m: model 
        @type  m: PDBModel
        @param f_pdb_out: name of pdb file to write
        @type  f_pdb_out: str
        """
        ## Whatif deletes hydrogens anyway
        m = m.compress( m.maskHeavy() )
        m.renumberResidues()
        m.writePdb( f_pdb_out, ter=0 )


    def cleanup( self ):
        """
        Tidy up the mess you created.
        """
        Executor.cleanup( self )

        if not self.debug:
            T.tryRemove( self.f_pdb )
            T.tryRemove( self.f_relativeASA )
            T.tryRemove( self.f_residueASA )

            T.tryRemove('FOR???.DAT', wildcard=1)
            T.tryRemove('pdbout.tex')
            T.tryRemove('pdbout.txt')
            T.tryRemove('TEXSTORE.DAT')
            T.tryRemove('TEXTABLE.DAT')


    def isFailed( self ):
        """
        Overrides Executor method
        """
        return self.error is None 


    def finish( self ):
        """
        Overrides Executor method
        """
        Executor.finish( self )
        self.result = self.parse_whatif( )


    def __read_residueASA( self ):
        """
        Read solvent accessibility calculated with WHATIF and return array
        of ASA-values. First column is the total ASA, second column the ASA
        of the backbone, the third column is the ASA of the side-chain.

        @return: array (3 x len(residues)) with ASA values
        @rtype: array
        """
        ## [14:-27] skip header and tail lines line

        lines = open(self.f_residueASA).readlines()[14:-27]
        ASA_values = map(lambda line: string.split(line)[-3:], lines)
        ASA_values = N.array(map(lambda row: map(string.atof, row),
                                 ASA_values))

        return ASA_values


    def __read_relativeASA( self ):
        """
        Read the relative ASA data. Return dic indexed by resNumber
        with atom and relative accessibility (%)

        @return: relative atom ASA dictionary. e.g::
                 { 1: { 'N':34.6, 'CA':81.5 .. }, 2 : {.. }.... }
        @rtype: dict of dict
        """
        AtomASA = {} 

        logFile = open( self.f_relativeASA  )

        while 1:
            line = logFile.readline()

            # first line of new residue
            if line[:8] == 'Residue:':
                # get residue number
                residue = int( string.split(line)[1] )
                AtomASA[residue] = {}

                # skip two lines
                line = logFile.readline() 
                line = logFile.readline() 
                # read first line of atom list
                line = logFile.readline()
                while line[:20] != '                    ':

                    atom = string.split(line)[0]
                    value = float( string.split(line)[-1] )
                    AtomASA[ residue ][ atom ] = value

                    line = logFile.readline()

            if line == '': #eof
                break

        return AtomASA


    def __exposedResidues( self, ASA_values, sidechainCut=0.0,
                         backboneCut=0.0, totalCut=0.0  ):
        """
        Decide what is a surface exposed residue and what is not.
        sidechainCut, backboneCut, totalCut - float, cutoff value
        for what will be considered as a exposed residue. All
        three values have to pass the test.

        @param ASA_values: array with ASA values for side chains, backbone
                           and total calculated in L{__read_residueASA}.
        @type  ASA_values: array
        @param sidechainCut: cutoff ASA value for considering the side chain
                             to consider thew residue being exposed
                             (default: 0.0) 
        @type  sidechainCut: float
        @param backboneCut: cutoffvalue for back bone ASA
        @type  backboneCut: float 
        @param totalCut: cutoff for total ASA
        @type  totalCut: float   

        @return: residue mask, where 0 = burried
        @rtype: [1|0]
        """
        col_0 = N.greater( N.transpose(ASA_values)[0], totalCut )
        col_1 = N.greater( N.transpose(ASA_values)[1], backboneCut )
        col_2 = N.greater( N.transpose(ASA_values)[2], sidechainCut )

        col_012 = N.concatenate( ([col_0],[col_1],[col_2]) ) 

        exposedList = N.greater(N.sum(col_012), 0)

        return exposedList


    def parse_whatif( self ):
        """
        Takes a PDBModel and returns three lists::

          relativeASA - a list of atom relative accessibilities
                        (in %  - compared to the resideu in a
                        Gly-XXX-Gly tripeptide, an approximation
                        of the unfolded state)

                      -> N.array( 1 x N_atoms )

          residueASA  - array 3 x N_residues of float
                        acc.surf.area per res: sidechain, backbone, total


          burriedResidues - a burried atom mask
                            note: cutoff values set above

                          ->  [0, 1, 0, 0, 1, 1, 1, 0 ...]

        @return: relative accessibility, residue ASA and burried atom mask
        @rtype: array, array, [1|0]
        """
        # read ASA per residue log
        resASA = self.__read_residueASA( )

        # and decide what is a burried residue
        burriedResidues = self.__exposedResidues( resASA )

        # read the relative assesibility per atom log
        atomRelASA = self.__read_relativeASA( )

        return self.__asaAtomList( atomRelASA ), resASA, burriedResidues


    def __asaAtomList( self, relativeASA ):
        """
        convert ASA dict to list of values with len == atom number
        
        @param relativeASA: see __read_relativeASA()
        @type  relativeASA: dict

        @return: relative accessibilities e.g. [ 98.0, 0.0, 12.0, .. ],
                 Numpy array 1 x N(atoms)
        @rtype:  [float]
        """
        result = []

        resIndex = 1
        for res in self.model.resList():

            for atom in res:

                if relativeASA.has_key( resIndex ) \
                       and relativeASA[ resIndex ].has_key( atom['name'] ):

                    result.append( relativeASA[ resIndex ][ atom['name' ]] )
                else:
                    result.append( 0.0 )

            resIndex += 1

        return result


#############
##  TESTING        
#############
import Biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test Whatif interface"""

    TAGS = [ BT.OLD, BT.EXE ]  ## whatif is not officially supported any longer

    def test_Whatif(self):
        """Whatif test"""

        from Biskit import PDBModel

        ## Loading PDB...

        f = T.testRoot()+"/com/1BGS.pdb"
        m = PDBModel(f)

        m = m.compress( m.maskProtein() )
        m = m.compress( m.maskHeavy() )

        ## Starting WhatIf
        x = WhatIf( m, debug=0, verbose=0 )

        ## Running
        atomAcc, resAcc, resMask = x.run()

        if self.local:
            ## check that model hasn't changed
            m_ref = PDBModel(f)
            m_ref = m.compress( m.maskProtein() )
            for k in m_ref.atoms.keys():
                if not N.all(m_ref[k] == m[k]):
                    print 'Not equal ', k
                else:
                    print 'Equal ', k

            ## display exposed residues in PyMol
            from Pymoler import Pymoler
            pm = Pymoler()
            model = pm.addPdb( m, '1' )
            pm.colorRes( '1', resAcc[:,0] )
            pm.show()

            print "\nResult for first 10 atoms/residues: "
            print '\nAccessability (A^2):\n', atomAcc[:10]
            print '\nResidue accessability (A^2)'
            print '[total, backbone, sidechain]:\n', resAcc[:10]
            print '\nExposed residue mask:\n',resMask[:10]
            print '\nTotal atom    accessability (A^2): %.2f'%sum(atomAcc) 
            print '      residue accessability (A^2): %.2f'%sum(resAcc)[0]

        self.assertAlmostEqual( N.sum(resAcc[:,0]), 2814.6903, 7 ) 


if __name__ == '__main__':

    BT.localTest()






