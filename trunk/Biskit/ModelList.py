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

## last $Author$
## last $Date$
## $Revision$

"""
organise, sort, and filter list of PDBModels
"""
from DictList import DictList
from PDBModel import PDBModel

class ModelList( DictList ):
    """
    A list of diverse PDBModels. Special support is given to access and use
    entries in the info dictionaries for filtering and sorting.
    Some care is taken to avoid the adding of non-PDBModel objects (but
    somehow it is, certainly, still possible).

    See DictList for a complete description of all methods!

    Comparison to Trajectory:

    In contrast to Trajectory, ModelList does not require structures
    with identical atom content. This makes it of course less efficient and
    less powerful, in many respects. On the other hand, models are kept as
    seperate PDBModel instances with individual residue/atom profiles and
    info dictionaries.
    """

    def __init__( self, lst=[], item_type=PDBModel ):
        """
        @param lst: list
        @type  lst: [ dict ]
        @param item_type: type, class of allowed items [ dict ]
        @type  item_type: [ dict ]
        
        @raise raise DictListError: if list contains non-item_type item.
        """
        DictList.__init__( self, lst, item_type )


    def version( self ):
        """
        Version of class.
        
        @return: class version number
        @rtype: str
        """
        return DictList.version( self ) + '; ModelList $Revision$'


    def getItemValue( self, item, key, default=None ):
        """
        Get a value from a given item (dictionary). Overrides DictList method.
        
        @param item: possible entry of this list
        @type  item: dict
        @param key: dictionary key
        @type  key: any
        @param default: return value if key is not found (default: None)
        @type  default: any

        @return: any
        @rtype: any
        """
        return item.info.get( key, default )


    def file2model( self, file_or_model ):
        """
        Load PDBModel from file if necessary.
        
        @param file_or_model: file name or existing model
        @type  file_or_model: str OR PDBModel
        
        @return: existing model or model loaded from file
        @rtype: PDBModel
        """
        if isinstance( file_or_model, PDBModel ):
            return file_or_model

        return PDBModel( file_or_model )


    def _processNewItem( self, v, i ):
        """
        Called before an item is added to the list. Override but call.
        
        @param i: index
        @type  i: int
        @param v: value
        @type  v: dict

        @return: value
        @rtype: dict        
        """
        v = self.file2model( v )

        return DictList._processNewItem( self, v, i )


##
## TEST
##
if __name__ == '__main__':

    import Biskit.tools as T
    import glob, random

    f_lst = glob.glob( T.testRoot() + '/lig_pc2_00/pdb/*pdb.gz' )

    f_lst = f_lst[:5]

    print "Loading PDBs..."
    l = ModelList( f_lst )

    for m in l:
        m.info['score'] = random.random()

    print l.valuesOf('score')

    p = l.plot( 'score' )

    p.show()
