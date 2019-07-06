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
Base class for reading in Trajectory objects from different input formats.

This module provides helper classes for :class:`Biskit.md.Trajectory`. In most
cases, it should not be necessary to use it directly.

.. seealso:: `biskit.md.Trajectory`, `biskit.core.TrajParserFactory`,
"""

class TrajParserError( Exception ):
    """Error while parsing a structure file or object"""
    pass

class TrajParser:
    """
    **Abstract** base class for parsers that generate Trajectory objects from
    different source formats.
    """
    
    ## short free text description of the supported format
    description = ''
    
    def __init__(self, verbose=False, rmwat=False, analyzeEach=False):
        """
        Args:
            verbose (bool): print loading progress to STDERR
            rmwat (bool): remove water and other solvent molecules on the fly
            analyzeEach (bool): compare each frame's atom content to reference
        """
        self.verbose = verbose
        self.rmwat = rmwat
        self.analyzeEach = analyzeEach

    @staticmethod
    def supports( source ):
        """
        Override!

        The method is static and can thus be called directly with the parser
        class rather than with an instance:

        >>> if TrajParser.supports('myfile.crd'):
        >>>     ...

        Returns:
            bool: True if the given source is supported by this parser
                  implementation
        """
        raise NotImplementedError('issupported() is not implemented.')


    def parse2new( self, source, ref ):
        """
        Create a new Trajectory from the source.
        
        Args:
            source (str): file name or other input object
            ref (str or PDBModel): reference structure instance or file

        Returns:
           Biskit.Trajectory: new Trajectory instance
        """
        raise NotImplementedError('not implemented.')
    

##############
# Empty test #
##############

import biskit.test as BT

class Test(BT.BiskitTest):
    pass
