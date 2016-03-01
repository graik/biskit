## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py

#!/usr/bin/env python
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
## $Revision$
## last $Date$
## last $Author$

"""
Analyze HEX docking result.
"""

import Biskit.tools as t
from Biskit import Trajectory, mathUtils,  molUtils


import numpy as N
import copy

try:
    import biggles
except:
    biggles = 0

class AnalyzeError( Exception ):
    pass

class Analyzer:

    def __init__( self, verbose=1, **options ):
        """
        @param verbose: verbosity level (default: 1)
        @type  verbose: 1|0        
        @param options: needs::
                         rec,lig - file name, receptor, ligand trajectories
                         ref     - file name, pickled reference complex
        @type  options: any

        @raise AnalyzeError: if atoms are not aligned
        """
        self.options = options

        if verbose: t.flushPrint("\nLoading...")
        self.t_lig = t.load( options['lig'] )
        self.t_rec = t.load( options['rec'] )
        self.com= t.load( options['ref'] )

        ## delete H from all players
        self.t_lig.removeAtoms( self.t_lig.getRef().maskH() )
        self.t_rec.removeAtoms( self.t_rec.getRef().maskH() )

        self.com.rec_model.remove( self.com.rec_model.maskH() )
        self.com.lig_model.remove( self.com.lig_model.maskH() )
        self.com.lig_transformed = None

        ## equalize atom content of free (trajectory) and bound models
        if verbose: t.flushPrint('\nCasting...')

        bnd_rec = self.com.rec()
        bnd_lig = self.com.lig_model

        self.t_rec.sortAtoms()
        self.t_lig.sortAtoms()

        bnd_rec = bnd_rec.sort()
        bnd_lig = bnd_lig.sort()

        #m_bnd_rec, m_t_rec = bnd_rec.equalAtoms( self.t_rec.getRef() )
        #m_bnd_lig, m_t_lig = bnd_lig.equalAtoms( self.t_lig.getRef() )

        i_bnd_rec, i_t_rec = bnd_rec.compareAtoms( self.t_rec.getRef() )
        i_bnd_lig, i_t_lig = bnd_lig.compareAtoms( self.t_lig.getRef() )

        #self.t_rec.removeAtoms( N.logical_not( m_t_rec ) )
        #self.t_lig.removeAtoms( N.logical_not( m_t_lig ) )

        self.t_rec = self.t_rec.takeAtoms( i_t_rec )
        self.t_lig = self.t_lig.takeAtoms( i_t_lig )

        #self.mask_free_lig = m_t_lig
        #self.mask_free_rec = m_t_rec
        self.i_free_lig = i_t_lig
        self.i_free_rec = i_t_rec

        #bnd_rec.remove( N.logical_not( m_bnd_rec ) )
        #bnd_lig.remove( N.logical_not( m_bnd_lig ) )
        bnd_rec = bnd_rec.take( i_bnd_rec )
        bnd_lig = bnd_lig.take( i_bnd_lig )

        #self.mask_bnd_rec = m_bnd_rec
        #self.mask_bnd_lig = m_bnd_lig
        self.i_bnd_rec = i_bnd_rec
        self.i_bnd_lig = i_bnd_lig

        ## put 'equalized' models back into ref complex
        self.com.rec_model = bnd_rec
        self.com.lig_model = bnd_lig
        self.com.lig_transformed = None

        ## check
        if not self.t_rec.getRef().equals( self.com.rec() )[1] or \
           not self.t_lig.getRef().equals( self.com.lig() )[1]:

            raise AnalyzeError('Atoms are not aligned.')

        ## native contact matrix
        self.contacts = self.com.atomContacts()

        self.hexContacts = None


    def __categorizeHexSurf(self, cutoff=0.1):
        """
        Compare complexes of list to native complex to see if
        their contact surfaces overlapp with the native complex.
        
        @param cutoff: fraction cutoff for defining a overlap (default: 0.1)
        @type  cutoff: float
        
        @return: list of len(self.hexContacts) overlapping with
                 native contact surface of lig and rec (0 - no overlap,
                 1 - rec OR lig overlapps, 2- rec AND lig overlapps)
        @rtype: [0|1|2]
        """
        result = [ self.com.fractionNativeSurface( c, self.contacts )
                   for c in self.hexContacts ]

        result = [ N.sum( N.greater( o, cutoff ) ) for o in result ]
        return result


    def setHexComplexes(self, com_lst, n=20 ):
        """
        add contact matrices of hex-generated (wrong) complexes for comparison
        
        @param com_lst: ComplexList with contacts calculated
        @type  com_lst: ComplexList
        """
        t.flushPrint('adding hex-generated complexes')

        self.hexContacts = []

#        f = lambda a: molUtils.elementType(a['element']) == 'p'

        i = 0
        while i < n:

            com = com_lst[i]
            t.flushPrint('#')

            try:
                if com['fractNatCont'] == 0.0:#com['fnac_10'] == 0.0:
                    com.rec_model.remove( com.rec().maskH() )
                    com.lig_model.remove( com.lig_model.maskH() )

                    com.rec_model = com.rec_model.sort()
                    com.lig_model = com.lig_model.sort()

                    com.rec_model.keep( self.i_free_rec )
                    com.lig_model.keep( self.i_free_lig )
                    com.lig_transformed = None

                    self.hexContacts += [ com.atomContacts() ]

                else:
                    n += 1
                i+= 1
            except:
                print t.lastError()

        self.hexSurfaces = self.__categorizeHexSurf( 0.2 )


    def random_contacts( self, contMat, n, maskRec=None, maskLig=None ):
        """
        Create randomized surface contact matrix with same number of
        contacts and same shape as given contact matrix.
        
        @param contMat: template contact matrix
        @type  contMat: matrix
        @param n: number of matrices to generate
        @type  n: int
        @param maskRec: surface masks (or something similar)
        @type  maskRec: [1|0]
        @param maskLig: surface masks (or something similar)
        @type  maskLig: [1|0]
        
        @return: list of [n] random contact matricies
        @rtype: [matrix]
        """
        a,b = N.shape( contMat )
        nContacts = N.sum( N.sum( contMat ))

        if not maskLig:
            r_size, l_size = N.shape( contMat )
            maskLig = N.ones( l_size )
            maskRec = N.ones( r_size )

        c_mask = N.ravel( N.outer( maskRec, maskLig ) )
        c_pos = N.nonzero( c_mask )[0]

        # get array with surface positions from complex
        cont = N.take( N.ravel(contMat), c_pos )
        length = len( cont )

        result = []

        for i in range( n ):
            # create random array
            ranCont = mathUtils.randomMask( nContacts,length )

            # blow up to size of original matrix
            r = N.zeros(a*b)
            N.put( r, c_pos, ranCont)

            result += [ N.reshape( r, (a,b) ) ]

        return result


    def __shuffleList(self, lst ):
        """
        shuffle order of lst

        @param lst: list to shuffle
        @type  lst: [any]
        
        @return: shuffeled list
        @rtype: [any]
        """
        pos = N.random.permutation( len( lst ))
        return N.take( lst, pos )


    def shuffledLists( self, n, lst, mask=None ):
        """
        shuffle order of a list n times, leaving masked(0) elements untouched

        @param n: number of times to shuffle the list
        @type  n: int
        @param lst: list to shuffle
        @type  lst: [any]
        @param mask: mask to be applied to lst
        @type  mask: [1|0]

        @return: list of shuffeled lists
        @rtype: [[any]]        
        """
        if not mask:
            mask = N.ones( len(lst)  )

        if type( lst ) == list:
            lst = N.array( lst )
        
        pos = N.nonzero( mask )[0]

        rand_pos = N.array( [ self.__shuffleList( pos ) for i in range(n) ] )

        result = []
        for p in rand_pos:

            r = copy.copy( lst )
            N.put( r, p, N.take( lst, pos ) )
            result += [r]

        return result


    def __plotSequence(self, list, **arg):
        """
        add curve to current plot, with x1..xn = 0..n and y1..yn = |list|

        @param list: list to plot
        @type  list: [any]        
        """
        if not biggles:
            raise ImportError, 'biggles module could not be imported.'        

        self.plot.add( biggles.Curve( range( len(list) ), list, **arg ) )

        if arg.has_key('label'):
            self.plot.lastY -= .05

            self.plot.add( biggles.PlotLabel( self.plot.lastX, self.plot.lastY,
                                              arg['label'], **arg ))

    def report(self):
        """
        override for actual plotting
        """
        pass


    def initPlot(self):
        """
        override for plot creation
        """
        if not biggles:
            raise ImportError, 'biggles module could not be imported.'
        
        self.page = biggles.FramedPlot()
        self.plot = self.page


    def plotAll( self ):
        """
        Show plot
        """
        self.initPlot()
        self.report()
        self.page.show()



#############
##  TESTING        
#############
import Biskit.test as BT
        
class Test(BT.BiskitTest):
    """Hmmer test"""

    
    def prepare(self):
        import tempfile
        self.f_out = tempfile.mktemp( '_test_rec.traj' )

    def cleanUp(self):
        t.tryRemove( self.f_out )

    def test_Analyzer( self):
        """Dock.Analyzer test """
        from Biskit import Trajectory
        from Biskit.Dock import ComplexList

        ## create a minimal 1-frame receptor trajectory from a pdb file
        self.t_rec = Trajectory( [t.testRoot()+'/rec/1A2P.pdb'],
                                 verbose=self.local)
        t.dump( self.t_rec, self.f_out )

        ## load a complex list
        cl = t.load( t.testRoot() + '/dock/hex/complexes.cl')

        self.a= Analyzer( rec = self.f_out,
                          lig = t.testRoot()+'/lig_pcr_00/traj.dat',
                          ref = t.testRoot()+'/com/ref.complex',
                          verbose = self.local)

        ## shuffle this list five times
        shuff_lst = self.a.shuffledLists( 5, range(8) )

        ## create two random contact matrices
        rand_mat = self.a.random_contacts( cl[0].atomContacts(), 2 )
        
        self.assertEqual( N.shape(rand_mat[1]), (1075, 876) ) 
    
        
if __name__ == '__main__':
   
    BT.localTest()
