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
## last $Author: graik $
## last $Date: 2007-06-22 17:34:10 +0200 (Fri, 22 Jun 2007) $
## $Revision: 459 $
"""
Fetch a PDBModel from the remote or a local NCBI PDB database.

@see L{PDBModel}
@see L{PDBParserFactory}
"""
import numpy as N
import urllib, re, tempfile, os

import Biskit.tools as T
import Biskit.settings as settings
import Biskit as B
from PDBParser import PDBParser, PDBParserError
from PDBParseModel import PDBParseModel


class PDBParseNCBI( PDBParseModel ):

    ex_resolution = re.compile(\
        'REMARK   2 RESOLUTION\. *([0-9\.]+|NOT APPLICABLE)' )

    ## resolution assigned to NMR structures
    NMR_RESOLUTION = 3.5

    @staticmethod
    def supports( source ):
        """
        The method is static and can thus be called directly with the parser
        class rather than with an instance::

        >>> if ParsePDBModel.supports( model ):
        >>>     ...

        @return: True if the given source is supported by this parser
        @rtype: bool
        """
        r = isinstance( source, str )
        return r and len(source) == 4

    @staticmethod
    def description():
        """
        The method is static and can thus be called directly with the parser
        class rather than with an instance::

        >>> if ParsePDBModel.description():
        >>>     ...

        @return: short free text description of the supported format
        @rtype: str
        """
        return 'fetch PDB entry from NCBI'

    def getLocalPDBHandle( self, id, db_path=settings.pdb_path ):
        """
        Get the coordinate file from a local pdb database.

        @param id: pdb code, 4 characters
        @type  id: str
        @param db_path: path to local pdb database
                        (default: L{settings.pdb_path})
        @type  db_path: str

        @return: the requested pdb file as a file handle
        @rtype: open file handle

        @raise PDBParserError: if couldn't find PDB file
        """
        id = str.lower( id )
        filenames = [os.path.join( db_path, '%s.pdb' % id),
                     db_path + '/pdb%s.ent' % id,
                     db_path + '/%s/pdb%s.ent.Z' %( id[1:3], id ) ]

        for f in filenames:
            if os.path.exists( f ):
                ## gzipped pdb file
                if f[-3:]=='.gz':
                    return gzip.open(f)
                ## the gzip module doesn't handle .Z files
                ## doesn't return open file handle 
                elif f[-2:]=='.Z':
                    p = subprocess.Popen( [ 'gunzip', '-c', f ],
                                          stdout=subprocess.PIPE )
                    return p.communicate()[0]
                ## uncompressed
                else:
                    return open(f)

        raise PDBParserError( "Couldn't find PDB file locally.")


    def getRemotePDBHandle( self, id, rcsb_url=settings.rcsb_url ):
        """
        Get the coordinate file remotely from the RCSB.

        @param id: pdb code, 4 characters
        @type  id: str
        @param rcsb_url: template url for pdb download
                         (default: L{settings.rcsb_url})
        @type  rcsb_url: str

        @return: the requested pdb file as a file handle
        @rtype: open file handle

        @raise PDBParserError: if couldn't retrieve PDB file
        """
        try:
            from Bio import File
        except:
            raise PDBParserError('Could not find Biopython - ' + \
                                 'remote fetching of PDBs is not supported.')


        handle = urllib.urlopen( rcsb_url% (id,id) )

        uhandle = File.UndoHandle(handle)

        if not uhandle.peekline():
            raise PDBParserError( "Couldn't retrieve ", rcsb_url )

        return uhandle


    def parsePdbFromHandle(self, handle, first_model_only=True ):
        """
        Parse PDB from file/socket or string handle into memory.

        @param handle: fresh open file/socket handle to PDB ressource or string
        @type  handle: open file-like object or str
        @param first_model_only: only take first of many NMR models [True]
        @type  first_model_only: bool

        @return: pdb file as list of strings, dictionary with resolution
        @rtype: [str], {'resolution':float }
        @raise PDBParserError: if passed in string is too short
        """
        lines = []
        res_match = None
        infos = {}

        if type( handle ) is str:
            if len(handle) < 5000:
                raise PDBParserError( "Couldn't extract PDB Info." )
            handle =  cStringIO.StringIO( handle )

## if handle.peekline()[:6] != 'TITLE':
##     raise PDBParserError, 'Ressource does not seem to be a PDB:\n%r' %\
##     handle.peekline()

        for l in handle:
            lines += [ l ]

            res_match = res_match or self.ex_resolution.search( l )

            if first_model_only and l[:6] == 'ENDMDL':
                break
        
        if len(lines) < 10 and '<div>' in lines[0]:
            raise PDBParserError, 'No PDB found with this ID.'

        if res_match:
            if res_match.groups()[0] == 'NOT APPLICABLE':
                infos['resolution'] = self.NMR_RESOLUTION
            else:
                infos['resolution'] = float( res_match.groups()[0] )
        else:
            raise PDBParserError, 'No resolution record found in PDB.'

        return lines, infos


    def fetchPDB( self, id ):

        try:
            h = self.getLocalPDBHandle( id )
        except:
            h = self.getRemotePDBHandle( id )

        fname = tempfile.mktemp( '.pdb', 'ncbiparser_' )

        lines, infos = self.parsePdbFromHandle( h, first_model_only=1 )

        ## close if it is a handle
        try: h.close()
        except:
            pass

        f = open( fname, 'w', 1 )
        f.writelines( lines )
        f.close()

        m = B.PDBModel( fname )
        m.disconnect()
        m.pdbCode = id
        m.info.update( infos )

        T.tryRemove( fname )

        return m

    def update( self, model, source, skipRes=None, updateMissing=0, force=0,
                headPatterns=[]):
        """
        Update empty or missing fields of model from the source.

        Profiles that are taken from the source are labeled 'changed'=0.
        The same holds for coordinates (xyzChanged=0).
        However, existing profiles or coordinates or fields remain untouched.

        @param model: existing model
        @type  model: PDBModel
        @param source: PDB code
        @type  source: str
        @param skipRes: list residue names that should not be parsed
        @type  skipRes: [ str ]
        @param updateMissing: check source for additional profiles [0] 
        @type  updateMissing: 1|0
        """
        try:
            if force or updateMissing or self.needsUpdate( model ):

                s = self.fetchPDB( source )

                super( PDBParseNCBI, self ).update(
                    model, s, skipRes=skipRes, updateMissing=updateMissing,
                    force=force )

        except Exception, why:
            raise PDBParserError, "Cannot fetch PDB from %s, "\
                  % str(source) + "Reason:\n" + str(why)

        ## override source set by PDBParseModel
        model.source = source


#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test"""

    def test_PDBParseNCBI( self ):
        """PDBParseNCBI test"""

        ## loading output file from X-plor
        if self.local:
            print 'Loading pdb file ..'

        self.p = PDBParseNCBI()
        self.m = self.p.parse2new( '1A2P')

        self.assert_( len(self.m) == 3042 )
    
    def test_PDBParseNCBI_fail(self):
        """PDBParseNCBI wrong ID test"""
        ## loading output file from X-plor
        if self.local:
            print 'Requesting non-existing ID ..'
        
        self.p = PDBParseNCBI()
        
        with self.assertRaises(PDBParserError):
            self.p.parse2new('liv5')
       
        
    
if __name__ == '__main__':

    BT.localTest()
