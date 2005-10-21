#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##
## $Revision$
## last $Date$
## last $Author$

import Biskit.tools as t
from Biskit import Trajectory, mathUtils,  molUtils

import Numeric as N
import RandomArray
import copy

import biggles

class AnalyzeError( Exception):
    pass

class Analyzer:

    def __init__( self, **options ):
        """
        rec,lig - file name, receptor, ligand trajectories
        com     - file name, pickled reference complex
        """
        self.options = options

        t.flushPrint("\nLoading...")
        self.t_lig = t.Load( options['lig'] )
        self.t_rec = t.Load( options['rec'] )
        self.com= t.Load( options['ref'] )

        ## delete H from all players
        self.t_lig.removeAtoms( self.t_lig.getRef().maskH() )
        self.t_rec.removeAtoms( self.t_rec.getRef().maskH() )

        self.com.rec_model.remove( self.com.rec_model.maskH() )
        self.com.lig_model.remove( self.com.lig_model.maskH() )
        self.com.lig_transformed = None

        ## equalize atom content of free (trajectory) and bound models
        t.flushPrint('\nCasting...')
        
        bnd_rec = self.com.rec()
        bnd_lig = self.com.lig_model

        self.t_rec.sortAtoms()
        self.t_lig.sortAtoms()

        bnd_rec = bnd_rec.sort()
        bnd_lig = bnd_lig.sort()

        m_bnd_rec, m_t_rec = bnd_rec.equalAtoms( self.t_rec.getRef() )
        m_bnd_lig, m_t_lig = bnd_lig.equalAtoms( self.t_lig.getRef() )

        self.t_rec.removeAtoms( N.logical_not( m_t_rec ) )
        self.t_lig.removeAtoms( N.logical_not( m_t_lig ) )

        self.mask_free_lig = m_t_lig
        self.mask_free_rec = m_t_rec

        bnd_rec.remove( N.logical_not( m_bnd_rec ) )
        bnd_lig.remove( N.logical_not( m_bnd_lig ) )

        self.mask_bnd_rec = m_bnd_rec
        self.mask_bnd_lig = m_bnd_lig

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
        -> list of len(self.hexContacts) int
           overlapping with native contact surface of lig, rec yes,no
        """
        result = [ self.com.fractionNativeSurface( c, self.contacts )
                   for c in self.hexContacts ]

        result = [ N.sum( N.greater( o, cutoff ) ) for o in result ]
        return result


    def setHexComplexes(self, com_lst, n=20 ):
        """
        add contact matrices of hex-generated (wrong) complexes for comparison
        com_lst - ComplexList with Contacts
        """
        t.flushPrint('adding hex-generated complexes')
        
        self.hexContacts = []

        f = lambda a: molUtils.elementType(a['element']) == 'p'

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
                    
                    com.rec_model.remove( N.logical_not(self.mask_free_rec) )
                    com.lig_model.remove( N.logical_not(self.mask_free_lig) )
                    com.lig_transformed = None

                    self.hexContacts += [ com.atomContacts() ]
                    
                else:
                    n += 1
                i+= 1
            except:
                print t.lastError()

        self.hexSurfaces = self.__categorizeHexSurf( 0.2 )


    def random_contacts( self, contMat, n, maskRec, maskLig ):
        """
        Create randomized surface contact matrix with same number of
        contacts and same shape as given contact matrix.
        contMat - template contact matrix
        n - number of matrices to generate
        maskRec, maskLig - surface masks (or something similar)
        -> list of [n] random contact matricies
        """
        a,b = N.shape( contMat )
        nContacts = N.sum( N.sum( contMat ))

        c_mask = N.ravel( N.outerproduct( maskRec, maskLig ) )
        c_pos = N.nonzero( c_mask )

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
        """shuffle order of lst"""
        pos = RandomArray.permutation( len( lst ))
        return N.take( lst, pos )


    def shuffledLists( self, n, lst, mask ):
        """
        shuffle order of a list n times, leaving masked(0) elements untouched
        """
        if mask == None:
            mask == N.ones( len(lst)  )

        if type( lst ) == list:
            lst = N.array( lst )

        pos = N.nonzero( mask )

        rand_pos = N.array( [ self.__shuffleList( pos ) for i in range(n) ] )

        result = []
        for p in rand_pos:

            r = copy.copy( lst )
            N.put( r, p, N.take( lst, pos ) )
            result += [r]

        return result


    def __plotSequence(self, list, **arg):
        """add curve to current plot, with x1..xn = 0..n and y1..yn = |list|"""
        
        self.plot.add( biggles.Curve( range( len(list) ), list, **arg ) )

        if arg.has_key('label'):
            self.plot.lastY -= .05

            self.plot.add( biggles.PlotLabel( plot.lastX, plot.lastY,
                                              arg['label'], **arg ))

    def report(self):
        """override for actual plotting"""
        pass

    def initPlot(self):
        """override for plot creation"""
        self.page = biggles.FramedPlot()
        self.plot = self.page

    def plotAll( self ):
        self.initPlot()
        self.report()
        self.page.show()




