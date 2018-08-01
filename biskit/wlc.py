import numpy as N
import math

from biskit import gnuplot as G
from biskit import mathUtils as MU


class WormLikeChain:
    """
    Calculate polymer end-to-end distances and probabilities based on 
    worm-like chain model.
    
    Key terms in this model are:
        * contour length -- the (maximum) length of the fully stretched polymer
        * persistance length -- measure of chain flexibility, the minimum
          length needed to close a circle; for amino acids typically presumed
          to be around 4.5 A although values in literature vary.
    
    Most methods are provided in 2 versions: 
    (1) accepting absolut contour length in A as parameter -- e.g. r()
    (2) accepting length in amino acids as parameter  -- e.g. raa()
    """
    
    def __init__ ( self, caa=3.8, p=4.5, t=298):
        """
        @param caa: contour length per residue / amino acid in Angstroem [3.8]
        @type  caa: float
        @param p: persistence length in Angstroem [4.5]
        @type  p: float
        @param t: temperature in K [298]
        @type  t: float
        """
        self.caa = caa
        self.p = p
        self.t = t

    def r( self, c ):
        """
        Mean WLC end-to-end distance for given contour length.
        @param c: contour length in A
        @type  c: float
        @return: end-to-end distance in A
        @rtype: float
        """
        P = self.p ## persistence length
        return N.sqrt( 2*P*c*(1 - P/c*(1 - N.e**(-c/P)) ) )
        
    def raa( self, naa, extra=0. ):
        """
        WLC end-to-end distance depending on number of residues / amino acids.

        @param naa - int, number of amino acids
        @param extra - float, additional countour length in A

        @return: end-to-end distance in A
        @rtype: float
        """
        c = naa * self.caa  + extra ## contour length
        return self.r( c )
    
    def pr( self, c, r ):
        """
        Probability for a given end-to-end distance given a contour length.
        
        See: Zhou(2004): Polymer Models of Protein Stability, Folding, and Interactions
        
        @param c: contour length in A
        @type  c: float
        @param r: end-to-end distance in A
        @type  r: float
        @return: probability that the chain adopts given distance
        @rtype : float
        """
        p = self.p
        r2 = r**2
        c2 = c**2
        p2 = p**2
        r4 = r**4
        c4 = c**4
        
        P = 4*N.pi*r2 * (3./(4*N.pi*p*c))**3/2 * N.exp( -3.*r2/(4*p*c))\
          *(1- 5.*p/(4*c) - 2.*r2/c2 + 33.*r4/(80.*p*c**3) + 79.*p2/(160*c2)\
            + 329.*r2*p/(120*c**3) - 6799*r4/(1600.*c4) + \
            3441.*r**6/(2800*p*c**5) - 1089.*r**8/(12800.*p2*c**6))
        
        return P

    
    def praa( self, naa, r ):
        """
        Probability for a given end-to-end distance given length in residues.
        
        See: Zhou(2004): Polymer Models of Protein Stability, Folding, and Interactions
        
        @param naa: length in (amino acid) residues
        @type  naa: int
        @param r: end-to-end distance in A
        @type  r: float

        @return: probability that the chain adopts given distance
        @rtype : float
        """
        c = naa * self.caa 
        return self.pr( c, r )
    
    
    def force( self, r, c):
        """
        Worm Like Chain basic equation.

        @param r: end to end distance in A
        @type  r: float
        @param c: polymer contour length in A
        @type c: float

        @return: force in pN
        @rtype: float
        """
        p = self.p
        T = self.t

        k = 0.138 #pN * A

        Ec = (k * T )/(4.0*p)

        return Ec*(((1.-(r/c))**-2) - 1+(4*r/c))

    def contourLength(self, x, force) :
        """
        Returns polymer length for a given start to end separation and force.
        Very accurate solver for large values of x.

        @param x: average start to end distance (A)
        @type  x: float
        @param force: entropic force in pN
        @type force: float

        @return: contour length in A
        @rtype : float
        """
        T = self.t; p=self.p
        F= force

        k = 0.138 #K Boltzmann in pN  * A

        Ec = (k * T) /(4.0*p)

        ## x^3+ax^2+bx+c =0
        a = -2.*x*(F+3.*Ec)/F
        b = x*x*(9.*Ec+F)/F
        c = -4.*x*x*x*Ec/F

        return MU.cubic(a,b,c);


    def contourLength2(self, x, force) :
        """
        See getContourLengthFromX for description (from numerical recipes).
        It seems that this third grade solver is more accurate for small values 
        of x.
        @param x: average end to end distance (A)
        @type  x: float
        @param force: entropic force in pN
        @type force: float

        @return: contour length in A
        @rtype : float
        """
        T = self.t; p=self.p; 
        F= force

        k = 0.138 #pN A
        Ec = (k * T) /(4.*p)

        ## x^3+ax^2+bx+c =0
        a = -2.*x*(F+3*Ec)/F
        b = x*x*(9.*Ec+F)/F
        c = -4.*x*x*x*Ec/F

        ## Q and R are always real
        Q= (a*a -3.*b)/9.
        R = (2.*a*a*a - 9.*a*b +27.*c) / 54.

        if R*R < Q*Q*Q : 
            ## It has 3 real roots
            theta = math.acos(R/math.sqrt(Q*Q*Q))
            sqrQ =-2.* math.sqrt(Q)
            athi = a / 3.
            L1 = sqrQ*math.cos(theta/3)-athi
            L2 = sqrQ*math.cos((theta+6.28)/3.)-athi
            L3 = sqrQ*math.cos((theta-6.28)/3.)-athi
            return L1 , L2 , L3

        else:
            A = -N.sign(R)*(math.fabs(R) + math.sqrt(R*R-Q*Q*Q))**(1./3.)

            if A == 0:
                B = 0
            else:
                B = Q/A

            return (A+B)-(a/3.)

##############
## Test
##############
import biskit.test as BT

class Test(BT.BiskitTest):

    def test_WLC_dist2contour(self):
        """WLC distance to contour length"""
        p = 4.;lc = 5280.; x = 0.75*5280 

        wlc = WormLikeChain(p=p)

        f = wlc.force(x, lc)

        lc2 =  wlc.contourLength( x, force=f )
        lc3 =  wlc.contourLength2( x, force=f )

        self.assertEqual(lc2[0] - lc <2.,True)
        self.assertEqual(lc3[1] - lc <2.,True)
    
    def test_WLC_contour2dist(self):
        """WLC contour length to distance"""
        import biskit.gnuplot as G

        wlc = WormLikeChain()
        rc = [ (n, wlc.raa(n)) for n in range( 1, 40) ]
        
        if self.local:
            G.plot(rc) # plot most likely distance vs. length of chain in aa
            
        self.assertAlmostEqual(rc[12][1], 20.10, places=1)
    
    def test_WLC_contour2prob(self):
        """WLC contour length to probability"""
        import biskit.gnuplot as G
        
        wlc = WormLikeChain()
        pc = [ (r, wlc.praa(14,r)) for r in range( 1, 100) ]
        
        if self.local:
            G.plot( pc ) # plot probability vs. distance for 14aa long chain

        dpc = dict(pc)
        best = N.argmax(list(dpc.values()))
        self.assertEqual(best, 15)


if __name__ == '__main__':

    BT.localTest()
