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
Parse a certain file / memory object into a PDBModel.

This module provides helper classes for :class:`Biskit.PDBModel`. In most
cases, it should not be necessary to use it directly.

.. seealso:: `biskit.PDBModel`, `biskit.core.PDBParserFactory`,
             `biskit.core.PDBParseFile`, `biskit.core.PDBParseModel`,
             `biskit.core.PDBParsePickle`
"""
from biskit import StdLog, ErrLog
import biskit as B

class PDBParserError( Exception ):
    """Error while parsing a structure file or object"""
    pass

class PDBParser(object):
    """
    **Abstract** base class for parsers that generate PDBModel objects from
    different source formats.
    """

    def __init__( self, log=None ):
        """
        :param log: Log for warnings [default log to STDERR]
        :type  log: biskit.LogFile.LogFile

        Override if needed. Call parent method in overriding class!
        """
        self.log = log or ErrLog()


    @staticmethod
    def supports( source ):
        """
        Override!

        The method is static and can thus be called directly with the parser
        class rather than with an instance::

        >>> if PDBParser.supports('myfile.pdb'):
        >>>     ...

        :return: True if the given source is supported by this parser
                 implementation
        :rtype: bool
        """
        raise NotImplementedError('issupported() is not implemented.')


    @staticmethod
    def description():
        """
        Override!

        The method is static and can thus be called directly with the parser
        class rather than with an instance::

        >>> if PDBParser.description('myfile.pdb'):
        >>>     ...

        :return: short free text description of the supported format
        :rtype: str
        """
        raise NotImplementedError('description() is not implemented.')


    def update( self, model, source, skipRes=None, updateMissing=0, force=0,
                headPatterns=[] ):
        """
        Update empty or missing fields of model from the source. The
        model will be connected to the source via model.source.

        Override!

        Note for implementations:
        * Profiles that are taken from the source should be labeled
        'changed'=0 (regardless of their status in the source).
        * The same holds for coordinates (xyzChanged=0).
        * However, profiles or coordinates or fields existing in the model
          must remain untouched.

        :param model: existing model
        :type  model: PDBModel
        :param source: source PDB file or pickled PDBModel or PDBModel object
        :type  source: str || file || PDBModel
        :param skipRes: list residue names that should not be parsed
        :type  skipRes: [ str ]
        :param updateMissing: check source for additional profiles that are not
                           yet known for the model [False]
        :type  updateMissing: 1 || 0
        :param headPatterns: [(putIntoKey, regex)] extract given REMARKS
        :type  headPatterns: [(str, str)]
        """
        raise NotImplementedError('update() is not implemented.')


    def needsUpdate( self, model ):
        """
        This function is called by update() to decide whether or not to open
        the source.
        Override to make it more restrictive.

        :param model: model
        :type  model: PDBModel
        :return: true, if the model needs to be updated from its source
        :rtype: bool
        """
        return (model.xyz is None \
                or model.atoms.hasNoneProfile()\
                or model.residues.hasNoneProfile() )
##                or None in model.atoms.profiles.values() \
##                or None in model.residues.profiles.values() )


    def parse2new( self, source, disconnect=False, skipRes=None ):
        """
        Create a new PDBModel from the source.

        :param source: source PDB file or pickled PDBModel or PDBModel object
        :type  source: str || file || PDBModel
        :param disconnect: do *not* associate new model with the source [False]
        :type  disconnect: bool
        :param skipRes: list residues that should not be parsed
        :type  skipRes: [ str ]

        :return: new model (always of type PDBModel, regardless of source type)
        :rtype: PDBModel
        """
        m = B.PDBModel()
        self.update( m, source, updateMissing=True, skipRes=skipRes)

        if disconnect: m.disconnect()

        return m

##############
# Empty test #
##############

import biskit.test as BT

class Test(BT.BiskitTest):
    pass
