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

from Biskit.tools import *
import re
import tempfile
import Biskit.settings as settings

from Biskit import AmberRstParser 

def _use():

    print """
Concatenate 2 amber crd/vel files.

amberConcat.py -i sim1.crd sim2.crd -o sim_merged.crd -n n_atoms [-b
               -p |int_precission| -rst |last_sim1.rst last_sim2.rst| ]

    -n     number of atoms (obligatory)
    -b     traj has box info (3 additional coordinates)
    -p     when looking for overlapping block, round coord. to p positions
    -rst   try repairing last frame of sim1, sim2 from given restart file(s)
    """
    sys.exit( 0 )

class ParseError( Exception ):
    pass

class Concat:
    """
    Concatenate 2 Amber-generated crd files.
    """

    def __init__( self, f1, f2, out, n_atoms, boxInfo=0, precission=None,
                  lastRst=None ):

        self.f1 = absfile( f1 )
        self.f2 = absfile( f2 )
        self.o = absfile( out )
        self.n = n_atoms
        self.box = boxInfo
        self.precission = precission
        self.lastRst = lastRst

        self.lines_per_frame = self.n * 3 / 10

        if self.n % 10 != 0:
            self.lines_per_frame += 1
        if self.box:
            self.lines_per_frame += 1


    def line2numbers( self, l ):
        """convert line from crd/vel file to list of float numbers"""

        xnumber = "-*\d+\.\d+"              # optionally negtive number
        xspace  = ' *'                      # one or more space char
        ex = re.compile('('+xspace+xnumber+')')

        match = ex.findall( l )

        r = [ round(float(strCrd),3) for strCrd in match ] 

        return r


    def fuzzyEquals( self, l1, l2, prec ):
        """compare two lines of numbers allowing slight discrepancies."""

        if prec == None:
            return l1 == l2

        n1 = self.line2numbers( l1 )
        n2 = self.line2numbers( l2 )

        for i in range(10):
            if not (round(float(n1[i]), prec) == round(float(n2[i]), prec)):
                return 0

        return 1


    def go( self ):

        f1 = open( self.f1 )
        f2 = open( self.f2 )
        fo = open( self.o, 'w' )

        f2.readline()   ## skip first empty line of second crd
        f2_firstLine = f2.readline()

        ## get first line of repair restart block, if available
        rst_firstLine = None
        p = None
        if self.lastRst:
            print "Reading rst file for repair: ", stripFilename(self.lastRst)
            p = AmberRstParser( self.lastRst )
            rst_firstLine = p.getFirstCrdLine()

        ## write first empty line
        fo.write( f1.readline() )

        i = self.lines_per_frame
        frame = -1

        print "\nReading %s..." % stripFilename( self.f1 ),

        for l in f1:

            if i == self.lines_per_frame:
                i = 0
                frame += 1

                if frame % 100 == 0 and frame != 0:
                    flushPrint( '#' )

                if self.fuzzyEquals( l, f2_firstLine, self.precission ):
                    break

                if rst_firstLine and self.fuzzyEquals( l, rst_firstLine,
                                                       self.precission):
                    print "Replacing frame %i by frame from restart file" %frame
                    p.writeCrd( fo, lastAtom=self.n )
                    if self.box:
                        fo.write('   0.000   0.000   0.000\n' )
                    break
                    
            fo.write( l )
            i += 1

        print "%i frames. %i lines per frame. Last had %i lines." % \
              (frame, self.lines_per_frame, i+1)
        
        frame = 0
        i = self.lines_per_frame

        print "Reading %s..." % stripFilename( self.f2 ),

        fo.write( f2_firstLine )
    
        for l in f2:

            if i == self.lines_per_frame:
                i = 0
                frame += 1

                if frame % 100 == 0:
                    flushPrint( '#' )

            fo.write( l )
            i += 1

        print "%i frames. %i lines per frame. Last had %i lines." % \
              (frame, self.lines_per_frame, i+1)


        f1.close()
        f2.close()
        fo.close()

def concatMany( fout, n_atoms, box, precission, rst, f1, *fn ):
    """
    Hack: Use Concater to concat more than 2 files.
    """
    current_out = tempfile.mktemp('_concat.crd', dir=settings.tempDirShared)
    rst.reverse()
    fin = f1
    for f in fn:

        if f == fn[-1]:
            current_out = fout
        else:
            current_out = tempfile.mktemp('_concat.crd', dir=settings.tempDirShared)

        r = None
        if len( rst ) > 0:
            r = rst.pop()

        c = Concat( fin, f, current_out, n_atoms, box, precission, lastRst=r )

        c.go()

        fin = current_out
        

if __name__ == '__main__':

    if len( sys.argv ) < 3:
        _use()

    o = cmdDict( {} )
    f1 = o['i'][0]
    fn = o['i'][1:]
    n_atoms = int( o['n'] )

    precission = None
    if 'p' in o:
        precission = int( o['p'] )

    rst = toList( o.get('rst', [] ) )
    
    concatMany( o['o'], n_atoms, o.has_key('b'), precission, rst, f1, *fn )

##     c = Concat( f1, f2, o['o'], n_atoms, o.has_key( 'b' ), precission )
##     c.go()
