##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2006 Raik Gruenberg & Johan Leckner
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
## last $Author$
## last $Date$
## $Revision$
"""
Parse a PDB file into a PDBModel.

@see L{PDBModel}
@see L{PDBParserFactory}
"""
import Scientific.IO.PDB as IO
import Numeric as N

import Biskit.tools as T
import Biskit as B
from Biskit.PDBParser import PDBParser, PDBParserError


class PDBParseFile( PDBParser ):

    @staticmethod
    def supports( source ):
        """
        The method is static and can thus be called directly with the parser
        class rather than with an instance::

        >>> if PDBParser.supports('myfile.pdb'):
        >>>     ...

        @return: True if the given source is supported by this parser
                 implementation
        @rtype: bool
        """
        return (type(source) is str or isinstance(source, B.LocalPath)) and \
            (source[-4:].upper() == '.PDB' or
             source[-7:].upper() == '.PDB.GZ')


    @staticmethod
    def description():
        """
        The method is static and can thus be called directly with the parser
        class rather than with an instance::

        >>> if PDBParser.description('myfile.pdb'):
        >>>     ...

        @return: short free text description of the supported format
        @rtype: str
        """
        return 'PDB file'


    def idFromName( self, fname ):
        """
        Extract PDB code from file name.
        @param fname: file name
        @type  fname: str
        @return: first 4 letters of filename if available
        @rtype: str
        """
        name = T.stripFilename( fname )

        if len( name ) > 3:
            return name[:4]

        return ''


    def update( self, model, source, skipRes=None, lookHarder=0):
        """
        Update empty or missing fields of model from the source. The
        model will be connected to the source via model.source.

        @param model: existing model
        @type  model: PDBModel
        @param source: source PDB file
        @type  source: str
        @param skipRes: list residue names that should not be parsed
        @type  skipRes: [ str ]
        @param lookHarder: ignored
        @type  lookHarder: 1|0

        @raise PDBParserError - if something is wrong with the source file
        """

        try:
            ## atoms and/or coordinates need to be updated from PDB
            if self.needsUpdate( model ):

                atoms, xyz = self.__collectAll( source, skipRes )

                model.atoms = model.atoms or atoms

                model.xyz = model.xyz or xyz

                model.__terAtoms = model._PDBModel__pdbTer()

                model.fileName = model.fileName or source

                model.pdbCode = model.pdbCode or \
                                self.idFromName( model.fileName)
                               
        except:
            msg = self.__xplorAtomIndicesTest( source ) or ' '
            raise PDBParserError('Cannot read ' + str(source) + ' as PDB\n'\
                           '\ERROR: ' + T.lastError() + msg)

        model.setSource( source )


    def __xplorAtomIndicesTest( self, source ):
        """
        In some cases the setup with parallell xplor trajectories
        run out of atom indices when writing the pdb files to disc.
        When this happens (usualy for the TIP3 waters in the later
        of the 10 parallell trajectories) the atom indices get
        replaced with ***** which will cause the parsing to fail.
        The error message recieved is quite cryptic - this function
        is here to give a more comprehensible message.
        
        @param source: file that failed to be parsed
        @type  source: str
        """
        import re
        f = open( source, 'r' )
        lines = f.readlines()
        f.close()

        for i in range( len(lines) ):
            if re.match( '^ATOM\s{2}\*{5}', lines[i]):
                msg = """
Line %i to %i of the file %s contains invalid atom indices!

In some cases the setup with parallell xplor trajectories run out of atom indices when writing the pdb files to disc. When this happens (usualy for the TIP3 waters in the later of the 10 parallell trajectories) the atom indices get replaced with ***** which will cause the parsing to fail.

REMEDY: run the script fixAtomIndices.py
""" % (i, len(lines), source)

                return msg

    def __firstLetter( self, aName ):
        """
        Return first letter in a string (e.g. atom mane)

        @param aName: atom name
        @type  aName: str

        @return: first letter (i.e. not a number) from a string.
        @rtype: letter
        """
        try:
            i = int( aName[0] )
            return self.__firstLetter( aName[1:] )
        except:
            return  aName[0]


    def __collectAll( self, fname, skipRes=None ):
        """
        Parse ATOM/HETATM lines from PDB. Collect coordinates plus
        dictionaries with the other pdb records of each atom.
        REMARK, HEADER, etc. lines are ignored.

        Some changes are made to the dictionary from PDBFile.readline()::
            - the 'position' entry (with the coordinates) is removed
            - leading and trailing spaces are removed from 'name' ..
            - .. but a 'name_original' entry keeps the old name with spaces
            - a 'type' entry is added. Its value is 'ATOM' or 'HETATM'
            - a 'after_ter' entry is added. Its value is 1, if atom is
              preceeded by a 'TER' line, otherwise 0
            - empty 'element' entries are filled with the first non-number
              letter from the atom 'name'

        @param fname: name of pdb file
        @type  fname: str
        @param skipRes: list with residue names that should be skipped
        @type  skipRes: list of str

        @return: tuple of list of dictionaries from PDBFile.readline()
                 and xyz array N x 3
        @rtype: ( list, array )
        """
        items = []
        xyz   = []

        f = IO.PDBFile( fname )

        try:
            line, i = ('',''), 0

            while line[0] <> 'END' and line[0] <> 'ENDMDL':

                i += 1
                try:
                    line = f.readLine()
                except ValueError, what:
                    self.log.add('Warning: Error parsing line %i of %s' %
                                 (i, T.stripFilename( fname )) )
                    self.log.add('\tError: '+str(what) )
                    continue

                ## preserve position of TER records
                newChain = line[0] == 'TER'
                if newChain:
                    line = f.readLine()

                if (line[0] in ['ATOM','HETATM'] ):

                    a = line[1]

                    if skipRes and a['residue_name'] in skipRes:
                        continue

                    a['name_original'] = a['name']
                    a['name'] = a['name'].strip()

                    a['type'] = line[0]
                    if newChain: a['after_ter'] = 1

                    if a['element'] == '':
                        a['element'] = self.__firstLetter( a['name'] )

                    if a['position'].is_vector:
                        lst = [ a['position'][0],
                                a['position'][1],
                                a['position'][2]]
                        xyz.append( lst )
                    else:
                        xyz.append( a['position'] )

                    del a['position']

                    items += [ a ]

        except:
            raise PDBParserError("Error parsing file "+fname+": " + T.lastError())

        try:
            f.close()
        except:
            pass

        if len( xyz ) == 0:
            raise PDBParserError("Error parsing file "+fname+": "+
                            "Couldn't find any atoms.")

        return items, N.array( xyz, N.Float32 )
    

class Test:
    """
    Test class
    """

    def run( self, local=0 ):
        """
        run function test
        
        @param local: transfer local variables to global and perform
                      other tasks only when run locally
        @type  local: 1|0
        
        @return: coordinates of center of mass
        @rtype:  array
        """

        ## loading output file from X-plor
        if local:
            print 'Loading pdb file ..'

        p = PDBParseFile()
        m = p.parse2new( T.testRoot()+'/rec/1A2P.pdb')

        if local:
            globals().update( locals() )

        return N.sum( m.centerOfMass() )


    def expected_result( self ):
        """
        Precalculated result to check for consistent performance.

        @return: coordinates of center of mass
        @rtype:  array
        """
        return N.sum( N.array([ 29.53385022,  46.39655482,  37.75218589]))
        

if __name__ == '__main__':

    test = Test()

    assert abs( test.run( local=1 ) - test.expected_result() ) < 1e-8
