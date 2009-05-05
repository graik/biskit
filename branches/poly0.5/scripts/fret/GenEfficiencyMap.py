#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2009 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 2 of the
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
## last $Author$
## $Date: 2009-03-08 19:10:27 +0100 (Sun, 08 Mar 2009) $
## $Revision$

# MAIN AUTHOR: Victor Gil Sepulveda
# DESCRIPTION: Shows the efficiency map for a pair of chromophores.One is fixed and the other gets all 
#			the vector values given by the center of an unit radius sphere and its surface points, spaced
#			5 degrees. Distance vector used is (0,d,0). Chromophore 1 orientation is (0,0,1).
#			It also saves the data values as  spherescriptdata and the gnuplot script used as plotspherescript
#			for posterior plotting.
#			For extra features please look at FRETUtils.py (create3DFRETEfficiencySphere,plot3DFRETEfficiencySphere)
# STATUS: Finished and tested
# TESTS: 	- No numeric extinction coefficient
#		- Final map plot
#		- File creation


from Biskit.fret.fretutils import create3DFRETEfficiencySphere,plot3DFRETEfficiencySphere,sphericalVectors
from Biskit.fret.fret import FRET
import sys
from string import atof

if __name__ == '__main__':
	
	# Check if arguments are ok
	if len(sys.argv)<5:
		print "\nUsage: exe arg1 arg2 arg3 arg4 arg5\n"
		print "arg1->  Integral overlap (nm^6)"
		print "arg2->  Donor quantum yield"
		print "arg3->  Acceptor maximum extinction coefficient (M^-1 cm^-1)"
		print "arg4->  Distance between chromophores (A)"
		print "arg5->	 Output png file (without extension)\n"
		print "Note->	 For extra features please look at FRETUtils.py (create3DFRETEfficiencySphere,plot3DFRETEfficiencySphere)"
		
		exit( -1 )
	
	
	try:
		overlap = atof(sys.argv[1])
	except:
		print "Error: Impossible to do string to float conversion for arg. 1 (overlap)"
		exit( -1 )
	
	try:
		qyD = atof(sys.argv[2])
	except:
		print "Error: Impossible to do string to float conversion for arg. 2 (donor quantum yield)"
		exit( -1 )
	
	try:
		e_cof = atof(sys.argv[3])
	except:
		print "Error: Impossible to do string to float conversion for arg. 3 (extinction coefficient)"
		exit( -1 )
		
	try:
		d = atof(sys.argv[4])
	except:
		print "Error: Impossible to do string to float conversion for arg. 3 (distance)"
		exit( -1 )
	
	MapName = sys.argv[5]
	
	f = FRET( overlap, qyD,e_cof)

	acceptors , cycle = sphericalVectors( )

	#~ print acceptors

	results  = create3DFRETEfficiencySphere(  f, [0.,d,0.] ,[0.,0.,1.], acceptors  )

	#~ print results

	plot3DFRETEfficiencySphere(results,cycle,MapName,True,("set view map\n",))
