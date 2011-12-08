import math
import mathUtils as MU
import numpy as N


class WormLikeChainModel:

    def __init__ ( self, lc=0., p=0., t=298, f=0.):
        """
        @param lc: contour length in Angstroem
        @type  lc: float
        @param p: persistence length in Angstroem
        @type  p: float
        @param t: temperature in K [298]
        @type  t: float
        @param f: force in pN [0.]
        @type  f: float
        """
        self.lc = lc
        self.p = p
        self.t = t
        self.f = f


    def contourLength(self, x) :
        """
        Returns polymer length for a given start to end separation.
        Very accurate solver for large values of x.

        @param x: average start to end distance (A)
        @type  x: float

        @return: contour length in A
        @rtype : float
        """
        T = self.t; p=self.p; F= self.f

        k = 0.138 #K Boltzmann in pN  * A

        Ec = (k * T) /(4.0*p)

        ## x^3+ax^2+bx+c =0
        a = -2.*x*(F+3.*Ec)/F
        b = x*x*(9.*Ec+F)/F
        c = -4.*x*x*x*Ec/F

        return MU.cubic(a,b,c);


    def contourLength2(self, x) :
        """
        See getContourLengthFromX for description (from numerical recipes).
        It seems that this third grade solver is more accurate for small values 
        of x.
        @param x: average start to end distance (A)
        @type  x: float

        @return: contour length in A
        @rtype : float
        """
        T = self.t; p=self.p; F= self.f

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



    def force( self, x):
        """
        Worm Like Chain basic equation.

        @param x: start to end distance (A)
        @type  x: float

        @return: Force (pN)
        @rtype: float
        """
        p = self.p
        Lc = self.lc
        T = self.t

        k = 0.138 #pN * A

        Ec = (k * T )/(4.0*p)

        return Ec*(((1.-(x/Lc))**-2) - 1+(4*x/Lc))


##############
## Test
##############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """ Test cases for Polyfret"""

    def prepare(self):
        pass

    def cleanUp( self ):
        pass

    def test_WLC(self):
        """WLC test cases"""
        p = 4.;Lc = 5280.; x = 0.75*5280 

        wlc = WormLikeChainModel(lc=Lc, p=p)

        F =  wlc.force(x)

        wlc.f = F

        Lc2 =  wlc.contourLength( x )
        Lc3 =  wlc.contourLength2( x )

        self.assertEqual(Lc2[0] - Lc <2.,True)
        self.assertEqual(Lc3[1] - Lc <2.,True)

if __name__ == '__main__':
    BT.localTest()


