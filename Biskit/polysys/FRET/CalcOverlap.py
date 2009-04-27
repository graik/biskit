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
## last $Author: graik $
## $Date: 2009-03-08 19:10:27 +0100 (Sun, 08 Mar 2009) $
## $Revision: 658 $

# MAIN AUTHOR: Victor Gil Sepulveda
# DESCRIPTION: Used for calculating the spectral integral overlap for 2 chromophores.
# STATUS: Finished and tested
# TESTS: 	- Non-existing data files.
#		- No numeric extinction coefficient
#		- Different sized wavelength range for data files
#		- Different range (starting point) for wavelenths
#		- Final value of calculation (It seems that python version has some round off errors so result
#			is less accurate)

from FRETUtils import overlapCalc
from numpy import array
import sys
from string import atof

if __name__ == '__main__':
	
	# Check if arguments are ok
	if len(sys.argv)<4:
		print "\nUsage: exe arg1 arg2 arg3\n"
		print "arg1->  Normalized (to the peak) donor emission spectra file"
		print "arg2->  Normalized (to the peak) acceptor absortion spectra file"
		print "arg3->  Acceptor maximum extinction coefficient\n"
		print "Note: Emission spectra file and absortion spectra file must to be defined for the same wavelength range."
		exit( -1 )
	
	#Open data files...
	try:
		donor_em = open(sys.argv[1],"r")
	except (IOError,EOFError):
		print "-Error: "+sys.argv[1]+" file doesn't exist."
		exit( -1 )
	
	try:
		acceptor_ab = open(sys.argv[2],"r")
	except (IOError,EOFError):
		print "-Error: "+sys.argv[2]+" file doesn't exist."
		exit( -1 )
	
	# ...read data...
	lines = donor_em.readlines()
	de = []
	wlde = []
	for l in lines:
		data = l.split()
		if len(data) != 2:
			print "-Error while reading emission data file."
			print "The format of every line must be: "
			print "Wavelength Emission_Value"
			exit( -1 )
		else:
			de.append(atof(data[1]))
			wlde.append(atof(data[0]))
	
	lines = acceptor_ab.readlines()
	ab = []
	wlab = []
	for l in lines:
		data = l.split()
		if len(data) != 2:
			print "-Error while reading absorbance data file."
			print "The format of every line must be: "
			print "Wavelength Absorbance_Value"
			exit( -1 )
		else:
			ab.append(atof(data[1]))
			wlab.append(atof(data[0]))
	
	
	# ... and close them
	donor_em.close()
	acceptor_ab.close()
	
	# Check if data is compatible with function call
	
	if len(wlde) != len(wlab):
		print "-Error: data files must be defined for the same wavelength range."
		exit( -1 )
	if wlde[0] != wlab[0]:
		print "-Error: data files must be defined for the same wavelength range."
		exit( -1 )
	
	try:
		e_cof = atof(sys.argv[3])
	except:
		print "Error: Impossible to do string to float conversion for arg. 3 (extinction coefficient)"
		exit( -1 )
	
	
	# Calculate overlap integral and show it 
	print "Overlap integral is: "+ str(overlapCalc(array(wlab),array(ab),array(de),e_cof))
	
	exit( 0 )