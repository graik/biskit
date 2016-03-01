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
## last $Author$
## last $Date$
## $Revision$
"""
Filter out templates before Modeller run.
"""

import Biskit as B
import numpy as N
import Biskit.mathUtils as M
from Biskit import StdLog, EHandler


class TemplateFilter(object):
    """
    Helper class for Biskit.Mod.Modeller
    """

    #: default z-value cutoff: templates that are more than z_cutoff standard deviations
    #: below the average similarity to the target sequence are disposed
    Z_CUTOFF = 0.2

    #: default sequence identity cutoff: templates that share less than
    #: ID_CUTOFF with the target sequence are disposed
    ID_CUTOFF = 0.4

    def __init__( self, aln_info, target_id='target',  log=None, verbose=0 ):
        """
        @param aln_info: CheckIdentities instance
        @type  aln_info: Biskit.Mod.CheckIdentities
        @param target_id: PIR id of the target sequence in the alignment
        @type  target_id: str
        """
        self.log = log or StdLog()
        self.verbose = verbose

        if not aln_info.result:
            aln_info.go()

        self.aln = aln_info.result

        self.target_id = target_id

        #: extract array of template identifiers
        self.templates = self.aln.keys()
        self.templates.remove( target_id )
        self.templates = N.array( self.templates )

        #: extract array of sequence identities to target
        self.identities = [ self.aln[k]['info_ID'][target_id]
                            for k in self.templates ]

        self.identities = N.array( self.identities )

        #: will hold mask for kicking out templates (1..keep template; 0..dispose)
        self.filter_mask = N.ones( len( self.templates ) )


    def get_filtered( self ):
        """
        @return: the identifiers of the templates that survived the filtering
        @rtype: [ str ]
        """
        if not N.sometrue(self.filter_mask):
            self.filter_mask[ N.argmax( self.identities ) ] = 1

        r = N.compress( self.filter_mask, self.templates ).tolist()

        if self.verbose:
            self.log.add( '%i templates were filtered out' % \
                          ( len(self.templates) - len(r)) )
        return r


    def filter_z( self, cutoff=None ):
        """
        Filter out templates that are further away from the target sequence
        than the average template.

        @param zcutoff: z-value cutoff (default: TemplateFilter.Z_CUTOFF)
        @type  zcutoff: float

        @return: a mask with 0 for every template that is zcutoff standard
                 deviations below the average similarity to the target
        @rtype: numpy.array
        """
        cutoff = cutoff or self.Z_CUTOFF

        avg = N.average( self.identities, 0 )
        sd  = M.SD( self.identities ) or 1e-10  ## replace 0 standard deviation
        z   = (self.identities - avg) / sd

        r = N.greater( z, -1. * cutoff )

        self.filter_mask = r * self.filter_mask

        if self.verbose:
            self.log.add('%i of %i templates fall through z-value filter.' %
                         (len( N.flatnonzero( r==0 ) ), len(self.templates) ) )

        return r


    def filter_id( self, cutoff=None):
        """
        Kick out all templates with less than cutoff sequence identity
        to the target. If no template is above the cutoff, the one with the
        best identity is retained.

        @param cutoff: identity cutoff
        @type  cutoff: float
        @return: a mask with 0 for every template to be kicked out
        @rtype: numpy.array
        """
        cutoff = cutoff or self.ID_CUTOFF

        r = N.greater( self.identities, cutoff )

        self.filter_mask *= r

## if max(r) == 0:
##     r[ N.argmax( self.identities ) ] = 1

        if self.verbose:
            self.log.add('%i of %i templates fall through identity filter.' %
                         (len( N.flatnonzero( r==0 ) ), len(self.templates) ) )
        return r


#########
# TEST
#########

import Biskit.test as BT

class Test(BT.BiskitTest):
    """
    Test class
    """

    def prepare(self):
        from CheckIdentities import CheckIdentities
        import Biskit.tools as T

        self.aln_info = CheckIdentities( alignment = T.testRoot() + \
                                         '/Mod/project/t_coffee/final.pir_aln' )

    def test_TemplateFilter(self):
        """Mod.TemplateFilter test"""
        self.f = TemplateFilter( self.aln_info, verbose=self.local )

        self.f.filter_z( cutoff=0.25 )
        self.f.filter_id()

        self.filtered = self.f.get_filtered()

        self.assertEqual( self.filtered, ['1K2H_A', '1M31_A', '1K8U_A'] )

if __name__ == '__main__':

    BT.localTest()
