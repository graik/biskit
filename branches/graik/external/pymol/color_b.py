# rlc color_b.py version 4.0

from pymol import cmd

def color_b(selection="all",min_b=None, max_b=None, ramp=0,rainbow=0,nbins=10):

	"""
	
AUTHOR 

	Robert L. Campbell

USAGE

	color_b(selection='sel',ramp=0 or 1, rainbow=0 or 1, nbins=10)

	This function allows coloring of a selection as a function
	of B-value, either in a range of hues from blue through magenta
	to red (rainbow=0, the default) or in the colors of the rainbow
	(rainbow=1).  The division of B-value ranges can either be as a
	histogram (equal-sized B-value increments leading to unequal numbers
	of atoms in each bin: ramp=0) or as a ramp of B-value ranges with an
	equal number of atoms in each group (ramp=1).

	usage:
	  from within PyMOL do "run color_b.py" to load the function definition.  
		Then you can use for example:

		    color_b('(c;a or c;b)',ramp=1,rainbow=1,nbins=30)

		to color chains A and B with the rainbow colors in 30 colors of equal 
		numbers of atoms in each color.

	"""

	ramp=int(ramp)
	rainbow=int(rainbow)
	nbins=float(nbins)
	if nbins == 1:
		print "\nWARNING: You specified nbins=1, which doesn't make sense...resetting nbins=10\n"
		nbins=10.
	print "RAMP, RAINBOW, NBINS,", ramp,rainbow, nbins

	m = cmd.get_model(selection)
	sel = []
	b_list = []

	if len(m.atom) == 0:
		print "No atoms selected"

	else:
		for i in range(len(m.atom)):
			b_list.append(m.atom[i].b)

			# create list of b-limits (bottom limit of bin)
			b_lim=[]

		if max_b == None:
			max_b = max(b_list)
		if min_b == None:
			min_b = min(b_list)
		print "Minimum and Maximum B-values: ", min_b, max_b

		if ramp:
			# color in bins of equal numbers of atoms
			b_list.sort()

			# subtract 0.1 from the lowest B in order to ensure that the single
			# atom with the lowest B value doesn't get omitted
			b_list[0] = b_list[0] - 0.1

			bin_num = int(len(b_list)/nbins)
			for j in range(nbins):
				sel.append(selection + " and b > " + str(b_list[j*bin_num]))
				#print "Color select: ",sel[j]

		else:
			# color in bins of equal B-value ranges
			# subtract 0.1 from the lowest B in order to ensure that the single
			# atom with the lowest B value doesn't get omitted
			min_b = min_b - 0.1
			bin_width = (max_b - min_b)/nbins
			for j in range(nbins):
				sel.append(selection + " and b > " + str(min_b + j*bin_width))
				#print "Color select: ",sel[j]

		if rainbow:
			col=[]
			coldesc=[]
			for j in range(nbins):
				# must append the str(sel[j]) to the color name so that it is unique 
				# for the selection
				coldesc.append('col' + str(sel[j]) + str(j))
				# create colors using hsv scale starting at blue(240) through red(0)
				# in intervals of 240/(nbins-1)  (the "nbins-1" ensures that the last
				# color is in fact red (0)
				#convert to rgb
				hsv = (240-j*240/(nbins-1),1.,1.)
				#print "HSV: ",hsv
				#convert to rgb and append to color list
				rgb = hsv_to_rgb(hsv)
				#print "RGB: ",rgb
				col.append(rgb)
				cmd.set_color("col" + str(sel[j]) + str(j),col[j])
				#print coldesc[j],col[j]

		else:
			col=[]
			coldesc=[]
			for j in range(nbins):
				coldesc.append('col' + str(sel[j]) + str(j))
				# create colors in a ramp from blue through magenta to red
				rgb = [min(1.0, j*2/nbins), 0.0, min(1.0, (nbins-j-1)*2/nbins)]
				col.append(rgb)
				cmd.set_color("col" + str(sel[j]) + str(j),col[j])
				#print coldesc[j],col[j]

		for j in range(nbins):
			#print "Color select: ",sel[j]
			cmd.color(coldesc[j],sel[j])

# allow calling without parentheses: color_hist_b <selection=>, <ramp= >,<rainbow= >,<nbins= >

cmd.extend("color_b",color_b)


def color_local(selection="all",min_b=0.1, max_b=1.1, ramp=0,rainbow=0,nbins=10):
	"""
	Color using a fixed max amd min b-values for color scaling.

	Estimated b-values for local gamma atom fluctuations:
	     min  b-factor: 0.1
	     mean b-factor: 0.35
	     max  b-factor. 1.1
	     Gly and Ala (that has no gamma atoms) have b-values of 0
	"""
	color_b(selection, min_b+0.1, max_b+0.1, ramp, rainbow, nbins)
	cmd.color('grey','b=0')
	cmd.extend("color_local",color_local)

def color_c_local(selection="all",min_b=0.001, max_b=0.7, ramp=0,rainbow=0,nbins=10):
	"""
	Color using a fixed max amd min b-values for color scaling.

	Estimated b-values for local gamma complimentarity fluctuations:
	     min  b-factor: 0
	     mean b-factor: ?
	     max  b-factor. 0.33
	     Gly and Ala (that has no gamma atoms) have b-values of 0
	"""
	color_b(selection, min_b+0.1, max_b+0.1, ramp, rainbow, nbins)
	cmd.color('grey','b=0.0')
	cmd.extend("color_c_local",color_c_local)
	
def hsv_to_rgb(hsv):
	# algorithm found at: http://www.cs.rit.edu/~ncs/color/t_convert.html

	h = float(hsv[0])
	s = float(hsv[1])
	v = float(hsv[2])

	if( s == 0 ) :
		#achromatic (grey)
		r = g = b = v

	else:
		# sector 0 to 5
		h = h/60.            
		i = int(h)
		f = h - i			# factorial part of h
		#print h,i,f
		p = v * ( 1 - s )
		q = v * ( 1 - s * f )
		t = v * ( 1 - s * ( 1 - f ) )

		if i == 0:
			(r,g,b) = (v,t,p)
		elif i == 1:
			(r,g,b) = (q,v,p)
		elif i == 2:
			(r,g,b) = (p,v,t)
		elif i == 3:
			(r,g,b) = (p,q,v)
		elif i == 4:
			(r,g,b) = (t,p,v)
		elif i == 5:
			(r,g,b) = (v,p,q)
		else:
			(r,g,b) = (v,v,v)
			print "error, i not equal 1-5"

	return [r,g,b]

