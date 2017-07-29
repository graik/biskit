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

"""
Cache of ExeConfig instances
"""

import threading
from Biskit import ExeConfig

class ExeConfigCache:
    """
    Keep track of previously loaded ExeConfig instances in order to avoid the
    multiple parsing of one and the same config file. The class contains
    only static methods -- it is not necessary to create an instance of it.

    Usage:
     >>> my_config = ExeConfigCache.get( 'my_program' )

    Should be thread-save, i.e. multiple threads can simultaniously access
    the cache (not tested).
    """

    LOCK = threading.Lock()
    CACHE= {}


    @staticmethod
    def get( name, reload=0, **kw ):
        """
        Get the ExeConfig instance for the given program.
        
        @param name: program name
        @type  name: str
        @param reload: force new instance (re-read configuration file)
                       (default: 0)
        @type  reload: 0|1
        @param kw: options for L{ Biskit.ExeConfig() }; no effect for
                   cached entries unless reload=1
        @type  kw: key=value
        
        @return: ExeConfig object
        @rtype: ExeConfig
        """

        ExeConfigCache.LOCK.acquire()

        if not name in ExeConfigCache.CACHE or reload:
            ExeConfigCache.CACHE[ name ] = ExeConfig( name, **kw )

        r = ExeConfigCache.CACHE[ name ]

        ExeConfigCache.LOCK.release()

        return r


    @staticmethod
    def reset():
        """
        Empty the cache.
        """

        ExeConfigCache.LOCK.acquire()

        ExeConfigCache.CACHE = {}

        ExeConfigCache.LOCK.release()


    @staticmethod
    def len():
        """
        Number of entries in the cache.
        
        @return: length of cache
        @rtype: int
        """
        return len( ExeConfigCache.CACHE )


#############
##  TESTING        
#############
import Biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test"""

    
    def test_ExeConfigCache( self):
        """ExeConfigCache test """

        self.x = ExeConfigCache.get( 'xplor' )
        
        if self.local:
            print self.x

        self.assert_( self.x.dat_found )
        

if __name__ == '__main__':

    BT.localTest()
