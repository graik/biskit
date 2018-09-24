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

"""
Cache of ExeConfig instances
"""

## allow relative imports when calling module by itself for testing (pep-0366)
if __name__ == "__main__" and __package__ is None:
    import biskit.exe; __package__ = "biskit.exe"

import threading
from biskit import BiskitError

from .exeConfig import ExeConfig

class ExeConfigCacheError(BiskitError):
    pass

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
    def get( name, reload=False, **kw ):
        """
        Get the ExeConfig instance for the given program.
        
        :param name: program name
        :type  name: str
        :param reload: force new instance (re-read configuration file)
                       (default: False)
        :type  reload: bool
        :param kw: options for :class:` Biskit.ExeConfig() `; no effect for
                   cached entries unless reload=True
        :type  kw: key=value
        
        The keyword arg `configpath=[str]` can be used to look for
        configuration files in other than the default biskit locations.
        
        :return: ExeConfig object
        :rtype: ExeConfig
        """

        if not ExeConfigCache.LOCK.acquire(timeout=15):
            raise ExeConfigCacheError('Failed to aquire singleton lock.')

        try:
            if reload or not name in ExeConfigCache.CACHE:
                ExeConfigCache.CACHE[ name ] = ExeConfig( name, **kw )
    
            r = ExeConfigCache.CACHE[ name ]

        finally:
            ExeConfigCache.LOCK.release()

        return r


    @staticmethod
    def reset():
        """
        Empty the cache.
        """
        if not ExeConfigCache.LOCK.acquire(timeout=15):
            raise ExeConfigCacheError('Failed to aquire singleton lock.')

        ExeConfigCache.CACHE = {}

        ExeConfigCache.LOCK.release()


    @staticmethod
    def len():
        """
        Number of entries in the cache.
        
        :return: length of cache
        :rtype: int
        """
        return len( ExeConfigCache.CACHE )


#############
##  TESTING        
#############
import biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test"""

    
    def test_ExeConfigCache( self):
        """ExeConfigCache test """
        import os.path

        self.x = ExeConfigCache.get( 'xplor' )
        
        if self.local:
            print(self.x)

        self.assertTrue( os.path.exists(self.x.dat) )
        

if __name__ == '__main__':

    BT.localTest()
