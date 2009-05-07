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
##
## last $Author$
## last $Date: 2009-03-09 09:10:46 +0100 (Mon, 09 Mar 2009) $
## $Revision$

# MAIN AUTHOR: Victor Gil Sepulveda
# STATUS: Work in Progress

"""
Support functions for FRET Module. 
"""

"""
Comments from Raik: Looks good!

* overlapCalc: no print statements in methods, please!
* plot3D...: should really use an Executor instance. os.system-style code
             is notoriously difficult to maintain.
* tests?
"""
## Example of Gnuplot Executor:

from Biskit import Executor

class Gnuplotter( Executor ):
	"""
	should work like this...
	See Xplorer.py for a full example
	
	You can provide the script as a string and the completed script will be
	piped into the program -- in this case you won't even touch the hard disk.
	"""

	def __init__( self, template='', **kw ):
		super( Gnuplotter, self).init( 'gnuplot', template=template, **kw)



def dbPre (str ,type='int',default=0,ehandler = None):
	"""
	Internal use function for database string conversion. 
	
	@param str: string number to be converted, read from db 
	@type  str: string
	@param type: 
	@type  type: { 'name' : list/array }
	@param infos:    dictionary of existing meta infos
	@type  infos:    { 'name' : { 'date' : ... } }
	"""
	
	if   'X' in str:
		return default

	if type == 'int':
		return int(str)
	elif type== 'float':
		return float(str)
	elif type == 'int_range':
		aux = str.split(',')
		return [int(aux[0]),int(aux[1])]
	else:
		if ehandler:
			ehandler.error(" Type not known ( "+type+" )." )
		return default


def overlapCalc ( wl,acc_spectra,don_spectra, e_cof ):
	"""
	Calculates the overlap integral for two spectra.
	
	@param wl: Iterable container with wavelength range (in nm).  
			Ex . [300,301,302]
	@type wl: float list
	@param acc_spectra: Iterable container with (normalized to the peak) acceptor absortion spectra. May
					have the same length that @wl.
					Ex. [0.3,1.0,0.4]
	@type acc_spectra: float list
	@param don_spectra: iterable container with (normalized to the peak) donor emision spectra. May
					have the same length that @wl.
					Ex. [1.0,0.5,0.4]
	@type don_spectra: float list
	@param e_cof: Extinction coefficient of the acceptor at maximum absorvance (in M^-1 cm^-1). 
				Ex. 43e3
	@type e_cof: float
	
	@return: Overlap integral 
	@rtype: float
	"""
	
	from numpy import sum
	
	area = sum(don_spectra)
	print  "Area: "+str(area)	
	don_spectra = don_spectra/area
	newintegral =sum(don_spectra)
	print "Now area is: "+str(newintegral)
	e_spectra =acc_spectra*e_cof
	return sum( wl**4*don_spectra*e_spectra)
	
def sphericalVectors( latdeg =0.087266462599716474  , longdeg = 0.087266462599716474):
	"""
	It returns a lists of vectors where each vector can be defined as the vector which starts in the center of an 
	unit radius sphere and a point of its surface.
	Surface points are generated every latdeg and longdeg degrees for latitude and longitude
	
	@param latdeg: Number of degrees for each latitude division. 
				Default value is 5 degrees, so defined sphere will have 360/5 latitude divisions.
	@type latdeg: float 
	@param longdeg: Number of degrees for each longitude division.
				Default value is 5 degrees, so defined sphere will have 360/5 longitude divisions.
	@type longdeg: float 
		
	@return: Overlap integral and lenght of a sphere contour in points
	@rtype: tuple {vector list, int}
	"""
	from emath import rotation
	from math import pi
	from numpy import array,matrix
	
	vectors = []
	others = []
	alpha = -pi/2.
	while alpha <pi/2.:
		beta = -pi/2.
		while beta < pi/2.:
			v = array([1.,0,0])*matrix(rotation(beta,alpha,0.))
			vectors.append((v[0,0], v[0,1],v[0,2]))
			others.append((-v[0,0], v[0,1],v[0,2]))
			beta+=latdeg
		alpha +=longdeg
		others.reverse()
		cycle = len(others)*2
		for v in others:
			vectors.append(v)
		others = []
	return vectors, cycle

def create3DFRETEfficiencySphere(  f = None,distance=(1,0,0),donor = (1,0,0), acceptors = []):
	"""
	Calculates the FRET efficiency of one chromophore trans. dipole moment with a list of other chromophore's 
	transition dipole moments
	
	@param f: FRET object for parameter storage and calculations
	@type latdeg: FRET 
	@param distance: Distance from the start of one chromophore dipole to the start of the dipole moment
				of the other chromophore.
	@type distance: float tuple
	@param donor: Transition dipole moment of one chromophore.
	@type donor: float tuple
	@param acceptors: List of transition dipole moments (of the other chromophore involved)
	@type acceptors: float tuple	
	
	@return: list of tuples containing the acceptor transition dipole moment (those in @acceptors) and its
		FRET Transfer Efficiency.
	@rtype: tuple list
	"""
	from emath import norm
	
	r = norm(distance)
	results = []
	
	for a in acceptors:
		f.calcK2( donor ,a, distance)
		results.append(( a[0],a[1],a[2],f.energyTransferEfficiency(r)))
	
	return results

def plot3DFRETEfficiencySphere( data = [] ,cycle = 0,filename= "",saveData = False ,more =("",),script = ("set hidden\n","set hidden3d\n","set pm3d\n","set pm3d depthorder\n","set ticslevel 0\n","set size square\n","splot \"spherescriptdata\" w l title \"Efficiency\"\n","pause 5\n",\
																							"set terminal png\n","set output \"myout\"\n","replot" )):
	"""
	Wrapper around gnuplot to make FRET Efficiency 3D plots. It will display it 5s by default, the save it to a file.
	
	@param data: list of tuples containing the acceptor transition dipole moment (those in @acceptors) and its
		FRET Transfer Efficiency.
	@type data: tuple list
	@param cycle:  Lenght of a sphere contour in points
	@type cycle: int 
	@param filename:  Name of PNG file for the plot to be saved
	@type filename: int 
	@param saveData:  If @True it also saves the gnuplot script and data used to do the plot ("plotspherescript" and
				"spherescriptdata").
	@type saveData: int 
	@param more:  First part of the gnuplot script.
	@type more: string iterable 
	@param script:  Gnuplot script to be used.
	@type script: string iterable
	"""
	
	import os
	
	# Gen. data file
	f = open("spherescriptdata",'w')
	c = 0
	for v in data:
		c +=1
		if c%cycle !=0:
			f.write(str(v[0])+" "+str(v[1])+" "+str(v[2])+" "+str(v[3])+"\n")
		else:
			f.write("\n")
	f.close()
	if not saveData:
		os.system("rm spherescriptdata")
	
	f = open("plotspherescript",'w')
	f.writelines(more + script)
	f.close()
	
	os.spawnlp(os.P_WAIT, 'gnuplot', 'gnuplot', 'plotspherescript')
	
	if filename != "":
		os.system("mv myout "+filename+".png")
	if not saveData:
		os.system("rm plotspherescript")
		os.system("rm myout")

