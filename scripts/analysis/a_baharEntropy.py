#!/usr/bin/env python
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
##
## Entropy estimate from free and bound XRay structures with the proteins
## treated as elastic lattice of CA.
## According to
## Bahar I, Atilgan AR, Erman B. (1997) Fold Des. 2(3):173-81
## the inverse of a special CA-contact matrix ( K-matrix ) is equivalent
## to the protein's fluctuation covariance matrix .
## Schlitter (1993) Chem. Phys. Lett. 215. 617
## uses such a covariance matrix (from MD) for an entropy estimate. 
##

from tools import *
from mathUtils import *
from PDBModel import *
import molUtils as MU

## Deprecated: would need to be converted to numpy -- not making any further changes
## without testing
from Numeric import *
from numpy.linalg import  svd, det
from LinearAlgebra import *

from gnuplot import *

class AnalyzeError( Exception ):
    pass

class BaharAnalyzer:
    """
    Take rec, lig, and com model of a complex. Calculate Entropy differences.
    """

    def __init__( self, **options ):

        self.options = options

        self.com = PDBModel( options['com'] )
        self.rec = PDBModel( options['rec'] )
        self.lig = PDBModel( options['lig'] )

        self.com = self.com.compress( self.com.maskCA() )
        self.rec = self.rec.compress( self.rec.maskCA() )
        self.lig = self.lig.compress( self.lig.maskCA() )

        ## chain number of rec and lig in com traj.
        self.cr = toIntArray( options.get('cr', 0) )
        self.cl = toIntArray( options.get('cl', 1) )


    def cast( self ):
        """Split com and cast free vs. bound"""

        self.com = self.com.sort()
        self.rec = self.rec.sort()
        self.lig = self.lig.sort()

        ## extract rec and lig coordinates from complex traj.
        rec_b = self.com.takeChains( self.cr )
        lig_b = self.com.takeChains( self.cl )

        if ( rec_b.lenAtoms() + lig_b.lenAtoms() != \
             self.com.lenAtoms() ):
            raise AnalyzeError('takeChain error.')

        mask_rec_f, mask_rec_b = self.rec.equalAtoms( rec_b )
        mask_lig_f, mask_lig_b = self.lig.equalAtoms( lig_b )

        self.rm_rec = len( nonzero(logical_not(mask_rec_f)) )
        self.rm_lig = len( nonzero(logical_not(mask_lig_f)) )

        ## apply to lig and rec
        self.rec = self.rec.compress( mask_rec_f )        
        self.lig = self.lig.compress( mask_lig_f )

        ## apply casting to complexTraj

        mask_com = concatenate( (mask_rec_b, mask_lig_b) )
        self.com = self.com.compress( mask_com )

        self.rm_com = len( nonzero(logical_not(mask_com)) )

        ## check, doesn't work if self.cr > self.cl
##         if not self.com.equals( self.rec.concat( self.lig ) )[1]:
##             raise AnalyzeError('Casting of complex failed.')


    def diagonalSum( self, m ):
        """replace diagonal of matrix m with column sum"""

        for i in range( len( m ) ):
            m[i,i] = 0

        s = sum( m )

        for i in range( len( m ) ):
            m[i,i] = - s[i]
        
        return m


    def kirchhofMatrix( self, m):
        """
        The inverse Kirchhof matrix is meant to be equivalent to the
        fluctuation covariance matrix as calculated from MD.
        m     - PDBModel
        -> array( m.lenAtoms() x m.lenAtoms() ) of int
        ij = 1, if d(ij) < 7 A and i != j
        ij = 0, else
        ij = sum( i ), if i==j
        """
        X = m.getXyz()

        r = pairwiseDistances( X, X )
        r = less( r, 7.0 ) * -1
        
        return self.diagonalSum( r )
    

    def getMatrices(self):
        """Calculate all necessary K-matrices"""

        ## kirchhof matrix for complex
        self.m_b = self.kirchhofMatrix( self.com )

        ## com matrix without rec-lig contacts
        rec_mask = [ x in self.cr for x in self.com.chainMap() ]
        lig_mask = [ x in self.cl for x in self.com.chainMap() ]

        rr_ll = outerproduct( rec_mask, rec_mask ) + \
                outerproduct( lig_mask, lig_mask )

        self.m_fcom = self.diagonalSum( self.m_b * rr_ll )

        ## matrices for free lig, rec
        self.m_rec = self.kirchhofMatrix( self.rec )
        self.m_lig = self.kirchhofMatrix( self.lig )


    def entropy( self, m ):
        """
        entropy( kirchhof_matrix) -> entropy [J/(K mol)]
        S is calculated as determinant = sum( ev ) of the co-variance matrix C.
        The K-matrix is the inverse of C, but with one zero eigenvalue.
        The determinant of K^-1 is hence assumed to be sum of its inverted
        eigenvalues, leaving out the zero one.
        """
##        e0 = svd( m)[1]
        e0 = eigenvectors(m)[0]

        ## remove near 0 eigenvalues
        e = compress( greater( e0, 1e-10 ), e0 )

        errWriteln("taking away %i zero eigenvalues." % (len(e0)-len(e)))
        
        return MU.boltzmann * MU.NA * 0.5 * sum( log( 1 / e ) )


    def fluctuation( self, m ):
        """
        The fluctuation is meant to be read from the diagonal of the inverse
        K-matrix.
        m - kirchhof matrix
        """
        e,v = eigenvectors( m )

        mask = greater( e, 1e-10 )
        e = compress( mask, e )
        v = compress( mask, v, 0)

        return sum( transpose(v)**2 / e, 1 )
    

    def entropyDiff( self, cast=1 ):
        """
        cast - remove non-equal atoms, (permanently alters internal models!)
        -> dS(bound - free), dS(bound - bound_no_contacts)
        """
        ent = self.entropy

        if cast:
            self.cast()

        self.getMatrices()

        ## bound - (free_rec + free_lig)
        ds_bf = ent( self.m_b ) - ent( self.m_rec ) - ent( self.m_lig )

        ## bound - bound_no_contacts
        ds_bb = ent( self.m_b ) - ent( self.m_fcom )

        return ds_bf, ds_bb


    def report( self ):
        """
        first 2 numbers are ent. diffs without casting, second 2 with.
        -> ( dS(bnd-free), dS(-bnd_no_cont), dS(bnd-free), sS(-bnd_no_cont))
        """
        return self.entropyDiff( 0 ) + self.entropyDiff( 1 )


    def reportLenghts( self ):
        """-> len_rec, len_lig, len_diff com - (rec+lig): (len in residues)"""
        lr, ll = self.rec.lenResidues(), self.lig.lenResidues()
        diff = self.com.lenResidues() - (lr + ll)

        return lr, ll, diff
        

###########################
## Analysis of many folders

def reportHeader():
    s = "\tlen AA,\t\tdiff b-f  castDiff \tdS before cast  dS after cast\n"
    s += "\t%3s\t%3s\t%4s\t%4s\t%4s\t%s\t%s\t%s\t%s\t%s" % \
           ('rec','lig','dcres','drec', 'dlig','dcom',
            'b-f', 'b-noco','b-f','b-noco')
    return s

def reportComplex( f ):

    errWriteln("Working on %s." % f )

    ## get -cr and -cl option for this complex
    o = file2dic( f+'/analysis_500-5.opt')

    o['com'] = f + '/com_wet/dry_com.model'
    o['rec'] = f + '/rec_wet/dry_rec.model'
    o['lig'] = f + '/lig_wet/dry_lig.model'

    a = BaharAnalyzer( **o )

    resLengths = a.reportLenghts()
    entropies = a.report()
    removed_atoms = (a.rm_rec, a.rm_lig, a.rm_com )

    return "%s\t%3i\t%3i\t%3i\t%2i\t%2i\t%2i\t%6.2f\t%6.2f\t%6.2f\t%6.2f" % \
           ((f,) + resLengths + removed_atoms + entropies )


def syntax():

    print """
a_baharEntropy.py -i |com_folder1 com_folder2 com_folder3..|

com_folder must contain com_wet/dry_com.model,
                        rec_wet/dry_rec.model,
                        lig_wet/dry_lig.model
                        and
                        analysis_500-5.opt for -cl and -cr option
E.g:
a_baharEntropy.py -i c11 c12 c13 > result.txt 2> log.txt
"""
    sys.exit(0)

if __name__ == '__main__':

    options = cmdDict( {} )

    if len( options ) < 1:
        syntax()
        
    print reportHeader()

    for f in toList( options['i'] ):

        try:
            print reportComplex( f )
        except:
            errWriteln("ERROR in %s." % f )
            errWriteln( lastErrorTrace() )
