##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2016 Raik Gruenberg & Johan Leckner
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
organise, sort, and filter list of PDBModels
"""
from DictList import DictList
from PDBModel import PDBModel
import Biskit.tools as T

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

        @raise DictListError: if list contains non-item_type item.
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


#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test class """

    def test_ModelList( self ):
        """ModelList test"""
        import random

        self.f_lst = [ T.testRoot() + '/rec/1A2P.pdb',
                       T.testRoot() + '/lig/1A19.pdb',
                       T.testRoot() + '/com/1BGS.pdb' ]

        ## Loading PDBs...
        self.l = ModelList( self.f_lst )

        self.seq = []
        for m in self.l:
            m.info['score'] = random.random()
            self.seq += [ m.compress( m.maskProtein() ).sequence() ]

        self.p = self.l.plot( 'score' )

        if self.local:
            print self.l.valuesOf('score')
            self.p.show()

        self.assertEqual( self.seq, self.EXPECTED )

    #: expected result
    EXPECTED = ['VINTFDGVADYLQTYHKLPDNYITKSEAQALGWVASKGNLADVAPGKSIGGDIFSNREGKLPGKSGRTWREADINYTSGFRNSDRILYSSDWLIYKTTDHYQTFTKIR', 'KKAVINGEQIRSISDLHQTLKKELALPEYYGENLDALWDCLTGWVEYPLVLEWRQFEQSKQLTENGAESVLQVFREAKAEGADITIILS', 'AQVINTFDGVADYLQTYHKLPDNYITKSEAQALGWVASKGNLADVAPGKSIGGDIFSNREGKLPGKSGRTWREADINYTSGFRNSDRILYSSDWLIYKTTDHYQTFTKIRKKAVINGEQIRSISDLHQTLKKELALPEYYGENLDALWDALTGWVEYPLVLEWRQFEQSKQLTENGAESVLQVFREAKAEGADITIILS']

if __name__ == '__main__':

    BT.localTest()
