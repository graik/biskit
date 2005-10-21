##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## organise, sort, and filter list of PDBModels

## last $Author$
## last $Date$
## $Revision$

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
        lst        - [ dict ]
        item_type  - type, class of allowed items [ dict ]
        !! raise DictListError, if list contains non-item_type item.
        """
        DictList.__init__( self, lst, item_type )

    def version( self ):
        return DictList.version( self ) + '; ModelList $Revision$'


    def getItemValue( self, item, key, default=None ):
        """
        Get a value from a given item (dictionary). Overrides DictList method.
        item    - dict, possible entry of this list
        key     - any, dictionary key
        default - any, return value if key is not found [None]
        -> any
        """
        return item.info.get( key, default )


    def file2model( self, file_or_model ):
        """
        Load PDBModel from file if necessary.
        file_or_model - str OR PDBModel, file name or existing model
        -> PDBModel, existing model or model loaded from file
        """
        if isinstance( file_or_model, PDBModel ):
            return file_or_model

        return PDBModel( file_or_model )

    def _processNewItem( self, v, i ):
        """
        Called before an item is added to the list. Override but call.
        i - int, index
        v - dict, value
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
