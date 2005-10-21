# rlc color_b.py version 4.0

from pymol import cmd

def color_q(selection="all",min_q=None, max_q=None, ramp=0,rainbow=0,nbins=10):

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
	  from within PyMOL do "run color_q.py" to load the function definition.  
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
	q_list = []

	if len(m.atom) == 0:
		print "No atoms selected"

	else:
		for i in range(len(m.atom)):
			q_list.append(m.atom[i].q)

			# create list of b-limits (bottom limit of bin)
			q_lim=[]

		if max_q == None:
			max_q = max(b_list)
		if min_q == None:
			min_q = min(b_list)
		print "Minimum and Maximum occupancy values: ", min_q, max_q

		if ramp:
			# color in bins of equal numbers of atoms
			q_list.sort()

			# subtract 0.1 from the lowest Q in order to ensure that the single
			# atom with the lowest Q value doesn't get omitted
			q_list[0] = q_list[0] - 0.1

			bin_num = int(len(q_list)/nbins)
			for j in range(nbins):
				sel.append(selection + " and q > " + str(q_list[j*bin_num]))
				#print "Color select: ",sel[j]

		else:
			# color in bins of equal Q-value ranges
			# subtract 0.1 from the lowest Q in order to ensure that the single
			# atom with the lowest q value doesn't get omitted
			min_q = min_q - 0.1
			bin_width = (max_q - min_q)/nbins
			for j in range(nbins):
				sel.append(selection + " and q > " + str(min_q + j*bin_width))
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
cmd.extend("color_q",color_q)


def color_global_pcr(selection="all",min_q=0.2, max_q=2.0, ramp=0,rainbow=0,nbins=10):
	"""
	Color using a fixed max and min b-values for color scaling.

	Estimated b-values for logal gamma atom fluctuations:
	     min  b-factor: 0.2
	     mean b-factor: 0.8
	     max  b-factor. 2.0
	"""
	color_q(selection, min_q+0.1, max_q+0.1, ramp, rainbow, nbins)
	cmd.extend("color_global_pcr",color_global_pcr)


def color_global_nopcr(selection="all",min_q=0.15, max_q=0.8, ramp=0,rainbow=0,nbins=10):
	"""
	Color using a fixed max and min b-values for color scaling.

	Estimated b-values for logal gamma atom fluctuations:
	     min  b-factor: 0.15
	     mean b-factor: 0.35
	     max  b-factor. 0.8
	"""
	color_q(selection, min_q+0.1, max_q+0.1, ramp, rainbow, nbins)
	cmd.extend("color_global_nopcr",color_global_nopcr)


def color_c_global_nopcr(selection="all",min_b=0.001, max_b=0.8, ramp=0,rainbow=0,nbins=10):
	"""
	Color using a fixed max amd min b-values for color scaling.

	Estimated b-values for global gamma complimentarity fluctuations:
	     min  b-factor: 0
	     mean b-factor: ?
	     max  b-factor. 0.4
	     Gly and Ala (that has no gamma atoms) have b-values of 0
	"""
	color_q(selection, min_b+0.1, max_b+0.1, ramp, rainbow, nbins)
	cmd.color('grey','q=0.0')
	cmd.extend("color_c_global_nopcr",color_c_global_nopcr)


def color_c_global_pcr(selection="all",min_b=0.001, max_b=1.8, ramp=0,rainbow=0,nbins=10):
	"""
	Color using a fixed max amd min b-values for color scaling.

	Estimated b-values for global gamma complimentarity fluctuations:
	     min  b-factor: 0
	     mean b-factor: ?
	     max  b-factor. 0.4
	     Gly and Ala (that has no gamma atoms) have b-values of 0
	"""
	color_q(selection, min_b+0.1, max_b+0.1, ramp, rainbow, nbins)
	cmd.color('grey','q=0.0')
	cmd.extend("color_c_global_pcr",color_c_global_pcr)


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

