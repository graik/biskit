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
"""Run WhatIf, calculate and collect accessibility data"""

import Numeric as N
import string
import tempfile

from Biskit import Executor, TemplateError
import Biskit.settings as S
import Biskit.tools as T


class WhatIf_Error( Exception ):
    pass

class WhatIf( Executor ):
    """
    Run WhatIf accessibility calculation
    
    pdbFileName - pdb file path
    relativeASA_log - whatif per atom relative ASA log file
    residueASA_log - whatif per residue absolute ASA log file
    whatifScriptName - whatif script file name

    Result: A WhatIf log file containing the accessibilities
            (in square Angstrom) per residue. The ASA is specified
            for the entire residue as well for the backbone
            and side chain.
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
        model - PDBModel

        ... and additional key=value parameters for Executor:
        debug    - 0|1, keep all temporary files                       [0]
        verbose  - 0|1, print progress messages to log     [log != STDOUT]
        node     - str, host for calculation (None->local) NOT TESTED [None]
        nice     - int, nice level                                     [0]
        log      - Biskit.LogFile, program log (None->STOUT)        [None]
        """
        Executor.__init__( self, 'whatif', template=self.whatif_script,
                           f_out='/dev/null', **kw )

        self.f_pdb = tempfile.mktemp('_whatif.pdb')
        self.f_relativeASA = tempfile.mktemp('_whatif_relative.log')
        self.f_residueASA = tempfile.mktemp('_whatif_residue.log')

        self.model = model.clone( deepcopy=1 )


    def prepare( self ):
        """Overrides Executor method."""
        self.__prepareModel( self.model, self.f_pdb ) 


    def __prepareModel( self, m, f_pdb_out ):
        """
        Don't change numbering of original model..
        """           
        ## Whatif deletes hydrogens anyway
        m = m.compress( m.maskHeavy() )
        m.renumberResidues()
        m.writePdb( f_pdb_out, ter=0 )


    def cleanup( self ):
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
        """Overrides Executor method"""
        return not self.error is None 


    def finish( self ):
        """Overrides Executor method"""
        Executor.finish( self )
        self.result = self.parse_whatif( self.output )
        
            
    def __read_residueASA( self ):
        """
        Read solvent accessibility calculated with WHATIF and return array
        of ASA-values. First column is the total ASA, second column the ASA
        of the backbone, the third column is the ASA of the side-chain.
        """
        ## [14:-27] skip header and tail lines line

        lines = open(self.f_residueASA).readlines()[14:-27]
        ASA_values = map(lambda line: string.split(line)[-3:], lines)
        ASA_values = N.array(map(lambda row: map(string.atof, row),
                                 ASA_values))

        return ASA_values


    def __read_relativeASA( self ):
        """
        Read the relative ASA data.
        Return dic indexed by resNumber with atom and
           relative accessibility (%)
           
        -> { 1: { 'N':34.6, 'CA':81.5 .. }, 2 : {.. }.... }
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

        -> [0, 1, 0, 0, 1, 1, 1, 0 ...] where 0 = burried
        """
        col_0 = N.greater( N.transpose(ASA_values)[0], totalCut )
        col_1 = N.greater( N.transpose(ASA_values)[1], backboneCut )
        col_2 = N.greater( N.transpose(ASA_values)[2], sidechainCut )

        col_012 = N.concatenate( ([col_0],[col_1],[col_2]) ) 

        exposedList = N.greater(N.sum(col_012), 0)

        return exposedList


    def parse_whatif( self, output ):
        """
        Takes a PDBModel and returns three lists:

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
        relativeASA - dict, see __read_relativeASA()
        
        -> [ 98.0, 0.0, 12.0, .. ], Numpy array 1 x N(atoms)
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


if __name__ == "__main__":
    
    from Biskit import PDBModel
    import Biskit.tools as T
    import glob

    print "Loading PDB..."

    f = glob.glob( T.testRoot()+'/lig_pc2_00/pdb/*_1_*pdb.gz' )[1]
    f = T.testRoot()+"/com_wet/1BGS.pdb"
    m = PDBModel(f)
    
    m = m.compress( m.maskProtein() )
    m.addChainFromSegid()
#    m = m.compress( m.chainMap() == 1 )
    m = m.compress( m.maskProtein() )
    m = m.compress( m.maskHeavy() )
    
    print "Starting WhatIf"
    x = WhatIf( m, debug=1, verbose=1 )

    print "Running"
    atomAcc, resAcc, resMask = x.run()

    print "\nResult for first 10 atoms/residues: "
    print '\nAccessability (A^2):\n', atomAcc[:10]
    print '\nResidue accessability (A^2)'
    print '[total, backbone, sidechain]:\n', resAcc[:10]
    print '\nExposed residue mask:\n',resMask[:10]
    print '\nTotal atom    accessability (A^2): %.2f'%sum(atomAcc) 
    print '      residue accessability (A^2): %.2f'%sum(resAcc)[0]

    m_ref = PDBModel(f)
    m_ref = m.compress( m.maskProtein() )
    for k in m_ref.atoms[0].keys():
        ref = [ m_ref.atoms[0][k] for i in range( m_ref.lenAtoms() ) ]
        mod = [ m.atoms[0][k] for i in range( m.lenAtoms() ) ]
        if not ref == mod:
            print 'Not equal ', k
        else:
            print 'Equal ', k
            
##     ## display exposed residues in PyMol
##     r = atomAcc
## #    mask = []
## #    for i in range(0, m.lenResidues() ):

## #        resAtoms = N.compress( N.equal( m.resMap(), i ), r )
## #        exposed = N.greater( resAtoms, 0.7)

## #        if N.sum(exposed) > len( resAtoms ) *1.0 / 4:
## #            mask += [1]
## #        else:
## #            mask += [0]

## #    r = m.residusMaximus(r)
## #    r = N.compress(m.maskCA(), r)
## #    r = N.greater(r, 40)
   
##     from Pymoler import *
##     pm = Pymoler()

##     model = pm.addPdb( m, '1' )
##     pm.colorRes( '1', r )

##     pm.show()
    

    



