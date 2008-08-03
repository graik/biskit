## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2007 Raik Gruenberg & Johan Leckner
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
## last $Author: graik $
## last $Date: 2007-06-22 17:34:10 +0200 (Fri, 22 Jun 2007) $
## $Revision: 459 $
"""
Parse a in-memory PDBModel instance into a new PDBModel

@see L{PDBModel}
@see L{PDBParserFactory}
"""
from Bio import File
import numpy.oldnumeric as N
import urllib, re, tempfile, os


import Biskit.tools as T
import Biskit.settings as settings
import Biskit as B
from PDBParser import PDBParser, PDBParserError
from PDBParseModel import PDBParseModel


class PDBParseNCBI( PDBParseModel ):

    ex_resolution = re.compile(\
            'REMARK   2 RESOLUTION\. *([0-9\.]+|NOT APPLICABLE)' )


    @staticmethod
    def supports( source ):
        """
        The method is static and can thus be called directly with the parser
        class rather than with an instance::

        >>> if ParsePDBModel.supports( model ):
        >>>     ...
        
        @return: True if the given source is supported by this parser
                 implementation (equivalent to isinstance( source, PDBModel) )
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
        return 'fetch PDB entry from NCBI web site or local database'

    def getLocalHandle( self, f ):
        """
        Open coordinate file from a local pdb database.

        @return: the requested pdb file as a file handle
        @rtype: open file handle

        @raise PDBParserError: if couldn't find PDB file
        """
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

    def getRemoteHandle( self, f ):
        """
        Open coordinate file via http socket.
        @return: the requested pdb file as url handle
        @rtype: open Bio.File.UndoHandle
        """

        handle = urllib.urlopen( f )
        uhandle = File.UndoHandle(handle)

        if not uhandle.peekline():
            raise PDBParseError( "Couldn't retrieve ", f )

        return uhandle


    def getRemotePDBAddress( self, id, rcsb_url=settings.rcsb_url ):
        """
        @param id: pdb code, 4 characters
        @type  id: str
        @param rcsb_url: template url for pdb download
                         (default: L{settings.rcsb_url})
        @type  rcsb_url: str
        """
        return rcsb_url% (id,id)
    

    def getLocalPDBAddress( self, id, db_path=settings.pdb_path ):
        """
        @param id: pdb code, 4 characters
        @type  id: str
        @param db_path: path to local pdb database
                        (default: L{settings.pdb_path})
        @type  db_path: str
        """
        id = id.lower()
        filenames = ['%s.pdb' % id,
                     db_path + '/pdb%s.ent' % id,
                     db_path + '/%s/pdb%s.ent.Z' %( id[1:3], id ) ]

        for f in filenames:
            if os.path.exists( f ):
                return f

        return None


    def parsePdbHeader(self, f ):
        """
	Parse PDB meta infos from open file/socket handle.

	@param handle: fresh open file/socket handle to PDB ressource or string
	@type  handle: open file-like object or str

        @return: pdb file as list of strings, dictionary with resolution
        @rtype: [str], {'resolution':float }
	@raise PDBParserError: if passed in string is too short
	"""
        res_match = None
        infos = {}

        if os.path.isfile( f ):
            handle = self.getLocalHandle( f )
        else:
            handle = self.getRemoteHandle( f )

        if type( handle ) is str:
            if len(handle) < 5000:
                raise PDBParserError( "Couldn't extract PDB Info." )
            handle =  cStringIO.StringIO( handle )

## 	if handle.peekline()[:6] != 'TITLE':
## 	    raise PDBParserError, 'Ressource does not seem to be a PDB:\n%r' %\
## 		  handle.peekline()

        for l in handle:

            res_match = res_match or self.ex_resolution.search( l )

            if res_match:
                if res_match.groups()[0] == 'NOT APPLICABLE':
                    infos['resolution'] = self.NMR_RESOLUTION
                else:
                    infos['resolution'] = float( res_match.groups()[0] )

                break
            
        if not res_match:
            raise PDBParserError, 'No resolution record found in PDB.'

        return infos


    def fetchPDB( self, id ):

        f = self.getLocalPDBAddress( id ) or self.getRemotePDBAddress( id )

        infos = self.parsePdbHeader( f )

        m = B.PDBModel( f )
        m.disconnect()

        return m


    def update( self, model, source, skipRes=None, updateMissing=0, force=0 ):
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
        model.setSource( source )


#############
##  TESTING        
#############
import Biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test"""

    def test_PDBParseModel( self ):
        """PDBParseNCBI test"""

        ## loading output file from X-plor
        if self.local:
            print 'Loading pdb file ..'

        self.p = PDBParseNCBI()
        self.m = self.p.parse2new( '1A2P')

        self.assert_( len(self.m) == 3042 )
                                

if __name__ == '__main__':

    BT.localTest()
