import numpy as N
import Biskit.gnuplot as G

class WormLikeChain:
    
    def __init__ ( self, caa=3.8, p=4.5, t=298, f=0.):
        """
        @param caa: contour length per residue / amino acid in Angstroem [3.8]
        @type  caa: float
        @param p: persistence length in Angstroem [4.5]
        @type  p: float
        @param t: temperature in K [298]
        @type  t: float
        @param f: force in pN [0.]
        @type  f: float
        """
        self.caa = caa
        self.p = p
        self.t = t
        self.f = f

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
        naa - int, number of amino acids
        extra - float, additional countour length in A
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
        c = naa * self.caa 
        return self.pr( c, r )
    
    def E( self, d ):
        """
        d - float, distance in A between fluorophores
        """
        R0 = 50
    
wlc = WormLikeChain()
rc = [ (n, wlc.raa(n)) for n in range( 1, 40) ]
G.plot( rc )

pc = [ (r, wlc.praa(14,r)) for r in range( 1, 100) ]
G.plot( pc )
