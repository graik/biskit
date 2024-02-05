## numpy-oldnumeric calls replaced by custom script; 09/06/2016
## Automatically adapted for numpy-oldnumeric Mar 26, 2007 by alter_code1.py

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2018 Raik Gruenberg & Johan Leckner
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
"""
Parse a PDB file into a PDBModel.

.. seealso:: `biskit.PDBModel` `biskit.PDBParserFactory`
"""
import biskit.core.scientificIO.PDB as IO
import biskit.core.oldnumeric as N0
import re

import biskit as B
import biskit.mathUtils as M
import biskit.tools as T
import biskit.biounit as BU
from biskit.core.pdbparser import PDBParser, PDBParserError
from biskit.core.localpath import LocalPath


class PDBParseFile( PDBParser ):

    #: default values for missing atom records
    DEFAULTS = {}
    
    #: default regular expressions for parsing REMARK entries
    RE_REMARKS = [ ('resolution', 
                    '2 RESOLUTION\. *([0-9\.]+|NOT APPLICABLE)' ),
                   
                   ]

    @staticmethod
    def supports( source ):
        """
        The method is static and can thus be called directly with the parser
        class rather than with an instance::

        >>> if PDBParser.supports('myfile.pdb'):
        >>>     ...

        :return: True if the given source is supported by this parser
                 implementation
        :rtype: bool
        """
        return (type(source) is str or isinstance(source, LocalPath)) and \
            (source[-4:].upper() == '.PDB' or
             source[-7:].upper() == '.PDB.GZ')


    @staticmethod
    def description():
        """
        The method is static and can thus be called directly with the parser
        class rather than with an instance::

        >>> if PDBParser.description('myfile.pdb'):
        >>>     ...

        :return: short free text description of the supported format
        :rtype: str
        """
        return 'PDB file'


    def idFromName( self, fname ):
        """
        Extract PDB code from file name.
        :param fname: file name
        :type  fname: str
        :return: first 4 letters of filename if available
        :rtype: str
        """
        name = T.stripFilename( fname )

        if len( name ) > 3:
            return name[:4]

        return ''


    def update( self, model, source, skipRes=None, updateMissing=0, force=0,
                headPatterns=[]):
        """
        Update empty or missing fields of model from the source. The
        model will be connected to the source via model.source.
        Profiles that are derived from the source are labeled 'changed'=0.
        The same holds for coordinates (xyzChanged=0).
        However, existing profiles or coordinates or fields remain untouched.

        :param model: existing model
        :type  model: PDBModel
        :param source: source PDB file
        :type  source: str
        :param skipRes: list residue names that should not be parsed
        :type  skipRes: [ str ]
        :param updateMissing: ignored
        :type  updateMissing: 1|0
        :param headPatterns: [(putIntoKey, regex)] extract given REMARKS
        :type  headPatterns: [(str, str)]

        :raise PDBParserError - if something is wrong with the source file
        """

        try:
            ## atoms and/or coordinates need to be updated from PDB
            if force or self.needsUpdate( model ):
    
                atoms, xyz, info = self.__collectAll( source, skipRes, 
                                                      headPatterns )
    
                keys = M.union( list(atoms.keys()),  list(self.DEFAULTS.keys()) )
    
                for k in keys:
    
                    a = model.atoms.get( k, default=0, update=False )
                    if (a is 0) or (a is None):
                    
                        dflt = self.DEFAULTS.get( k, None )
                        model.atoms.set(k, atoms.get(k, dflt), changed=0 )
    
                if model.xyz is None:
                    model.xyz = xyz
                    model.xyzChanged = 0
    
                model._resIndex  =None
                model._chainIndex=None
    
                model.fileName = model.fileName or source
    
                model.pdbCode = model.pdbCode or info.get('pdb_code', None) or \
                                self.idFromName( model.fileName)
                
                ## ## make biounit from the dictionary we have parsed                              if 'BIOMT' in info:
                ##     biomt = info['BIOMT']
                ##     model.biounit = BU.BioUnit(model, biomt)
                ##     del info['BIOMT']
                
                model.info.update( info )
                           
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
        
        :param source: file that failed to be parsed
        :type  source: str
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

        :param aName: atom name
        :type  aName: str

        :return: first letter (i.e. not a number) from a string.
        :rtype: letter
        """
        try:
            i = int( aName[0] )
            return self.__firstLetter( aName[1:] )
        except:
            return  aName[0]


    def __parseHeader( self, line_record ):
        """
        """
        return line_record[1]

    def __parseRemark( self, line_record, headPatterns=[]):
        """
        """
        l = line_record[1]
        for i in range( len(headPatterns)):
            key, regex = headPatterns[i]

            match = regex.search( l )
            if match:
                try:
                    value = None
                    del headPatterns[i]
                    value = match.groups()[0]
                    return { key : float( value ) }
                except:
                    return { key : value or l }

        return {}
        
    def __parseBiomt( self, pdbFile, firstLine):
        """
        Extract BIOMT (biological unit) information from REMARK 350 lines
        Creates a 'BIOMT' dictionary.
        """
        line = firstLine
        biomtDict = {}
        moleculeNum = -1

        while line[0] == 'REMARK' and line[1].startswith(' 350'):
            # 5 = len(' 350 ')
            biomtLine = line[1][5:].lstrip()
            if biomtLine.startswith('BIOMOLECULE:'): # start a new molecule

                if moleculeNum != -1:   
                    # lets update the dictionary with what we've got
                    biomtDict[moleculeNum] = (targetChains,rtList)

                #12 = len('BIOMOLECULE:')
                moleculeNum = int(biomtLine[12:].strip())
                targetChains = []
                rotation = []
                translation = []
                rtList = []

                matrixLine = 0

            if biomtLine.startswith('APPLY THE FOLLOWING TO CHAINS:'):  
            # parse targeted chains, we assume this comes after BIOMOLECULE line
                # 30 = len('APPLY THE FOLLOWING TO CHAINS:')
                targetChains.extend(c.strip() for c in biomtLine[30:].split(','))
            if biomtLine.startswith('AND CHAINS:'):  
                # 11 = len('AND CHAINS:')
                targetChains.extend(c.strip() for c in biomtLine[11:].split(','))

            if biomtLine.startswith('BIOMT'):  
            # parse rotate-translate matri{x/ces}, we assume this comes after BIOMOLECULE line
                matrixLine += 1
                # 6 = len('BIOMT#')
                rawCoords = biomtLine[6:].split()
                rotation.append([float(x) for x in rawCoords[1:4]])
                translation.append(float(rawCoords[4]))
                if matrixLine % 3 == 0:
                    rotation = N0.array( rotation )
                    translation = N0.transpose( [ translation ] )
                    rotation = N0.concatenate( (rotation, translation), axis=1 )
                    rtList.append(N0.array(rotation))
                    ## rtList.append((rotation,translation))
                    rotation = []
                    translation = []

            try:
                line = pdbFile.readLine()
            except ValueError as what:
                self.log.add('Warning: Error parsing line of %s' % T.stripFilename( str(pdbFile) ) )
                self.log.add('\tError: '+str(what) )
                continue
        # process last molecule group
        biomtDict[moleculeNum] = (targetChains,rtList)
        # return (indexed transformation dictionary , last line which isn't ours)
        return {'BIOMT': biomtDict}, line

    def __collectAll( self, fname, skipRes=None, headPatterns=[] ):
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

        :param fname: name of pdb file
        :type  fname: str
        :param skipRes: list with residue names that should be skipped
        :type  skipRes: list of str

        :return: tuple of (1) dictionary of profiles
                 and (2) xyz array N x 3
        :rtype: ( list, array )
        """
        xyz   = []

        aProfs = {}

        info = {}

        in_header = True
        
        headPatterns = headPatterns or self.RE_REMARKS
        patterns = [ (key, re.compile(ex)) for key,ex in headPatterns ]
        
        for k in B.PDBModel.PDB_KEYS:
            aProfs[k] = list()

        f = IO.PDBFile( fname )

        skipLine = False

        try:
            line, i = ('',''), 0

            while line[0] != 'END' and line[0] != 'ENDMDL':

                i += 1
                if not skipLine:
                    try:
                        line = f.readLine()
                    except ValueError as what:
                        self.log.add('Warning: Error parsing line %i of %s' %
                                     (i, T.stripFilename( fname )) )
                        self.log.add('\tError: '+str(what) )
                        continue
                else:
                    skipLine = False

                ## header handling
                if in_header and line[0] == 'HEADER':
                    info.update( self.__parseHeader( line ) )

                if in_header and line[0] == 'REMARK':
                    if line[1].startswith(' 350'):
                        biomtDict, line = self.__parseBiomt( f, line )
                        info.update( biomtDict )
                        # we've hogged a line beyond REMARK 350 records in 
                        # __parseBiomt(), now we need to process it here
                        skipLine = True
                        continue
                    else:
                        info.update( self.__parseRemark( line, patterns ) )
                    

                ## preserve position of TER records
                newChain = line[0] == 'TER'
                if newChain:
                    line = f.readLine()

                if (line[0] in ['ATOM','HETATM'] ):

                    if in_header: in_header = False  ## switch off HEADER parsing
                    
                    a = line[1]

                    if skipRes and a['residue_name'] in skipRes:
                        continue

                    a['name_original'] = a['name']
                    a['name'] = a['name'].strip()

                    a['type'] = line[0]

                    if newChain:
                        a['after_ter'] = 1
                    else:
                        a['after_ter'] = 0

                    if a['element'] == '':
                        a['element'] = self.__firstLetter( a['name'] )

                    xyz.append( a['position'] )

                    del a['position']

                    for k, v in a.items():
                        aProfs[k].append( v )

        except:
            raise PDBParserError("Error parsing file "+fname+": " \
                                 + T.lastError())
        try:
            f.close()
        except:
            pass

        if len( xyz ) == 0:
            raise PDBParserError("Error parsing file "+fname+": "+
                            "Couldn't find any atoms.")

        return aProfs, N0.array( xyz, N0.Float32 ), info
    
#############
##  TESTING        
#############
import biskit.test as BT
import time

def clock( s, ns=globals() ):
    import cProfile

    locals().update( ns )
    
    cProfile.run( s, 'report.out' )

    ## Analyzing
    import pstats
    p = pstats.Stats('report.out')
    p.strip_dirs()

    ## long steps and methods calling them
    p.sort_stats('cumulative').print_stats(20)
    p.print_callers( 20 )


class Test(BT.BiskitTest):
    """Test case"""

    def test_PDBParseFile( self ):
        """PDBParseFile test"""

        ## loading output file from X-plor
        if self.local:
            print('Loading pdb file ..')

        self.p = PDBParseFile()
        self.m = self.p.parse2new( T.testRoot()+'biounit/2V4E.pdb')
        if self.local:
            print(self.m.info)
        self.m.report( prnt=self.local,
                                plot=(self.local or self.VERBOSITY > 2) )
        self.m.biomodel(1)
        ##      self.m = self.p.parse2new( T.testRoot()+'/rec/1A2P_rec_original.pdb')
        ##      self.m2= self.p.parse2new( T.testRoot()+'/com/1BGS.pdb' )

        self.assertAlmostEqual( N0.sum( self.m.centerOfMass() ), 
                                -74.1017, 1 )



if __name__ == '__main__':

    BT.localTest()
