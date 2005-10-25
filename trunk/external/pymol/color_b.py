# rlc color_b.py version 6.2

import colorsys,sys
from pymol import cmd

# main function called from within PyMOL
def color_b(selection='all',item='b',mode='hist',gradient='bgr',nbins=11,sat=1.,value=1.,minimum='',maximum=''):
  """
    
  AUTHOR 

    Robert L. Campbell with enhancements from James Stroud

  USAGE

    color_b(selection='sel',item='b' or 'q'
      gradient='bgr' or 'rgb' or 'bwr' or 'rwb' or 'bmr' or 'rmb' or
      'rw' or 'wr' or 'gw' or 'wg' or 'bw' or wb' or 'gy' or 'yg' or 
      'gray' or 'reversegray'
      mode='hist' or 'ramp', [minimum=''],[maximum=20.],
      nbins=11, sat=1.0, value=1.0)

      The "item" argument allows specifying 'b' or 'q' as the item to
      color on.  The "color_q" function is the same as "color_b item=q".

      This function allows coloring of a selection as a function of
      B-value or occupancy, following a gradient of colours.  The
      gradients can be:

      'bgr': blue -> green   -> red
      'rgb': red  -> green   -> blue
      'bwr': blue -> white   -> red
      'rwb': red  -> white   -> blue
      'bmr': blue -> magenta -> red
      'rmb': red  -> magenta -> blue
      'rw' : red -> white
      'wr' : white -> red
      'gw' : green -> white
      'wg' : white -> green
      'bw' : blue -> white
      'wb' : white -> blue
      'gy' : green -> yellow
      'yg' : yellow -> green
      'gray' : black -> white
      'reversegray' : white -> black

      ('rainbow' and 'reverserainbow' can be used as synonyms for 
      'bgr' and 'rgb' respectively and 'grey' can be used as a synonym for 'gray').

      The division of B-value ranges can in either of two modes: 'hist' or
      'ramp'. 'hist' is like a histogram (equal-sized B-value increments
      leading to unequal numbers of atoms in each bin). 'ramp' as a ramp
      of B-value ranges with the ranges chosen to provide an equal number
      of atoms in each group.

      You can also specify the lower or upper limits of the data used to determine
      the color bins (minimum,maximum). e.g. color_b my_molecule, minimum=15., maximum=25.

      You can also specify the saturation and value (i.e. the "s" and "v"
      in the "HSV" color scheme) to be used for the gradient. The defaults
      are 1.0 for both "sat" and "value".

      In the case of the gray scale gradients, "sat" sets the minimum intensity 
      (normally black) and "value" sets the maximum (normally white)

    usage:
      from within PyMOL do "run color_b.py" to load the function definition.  
      Then you can use for example:

          color_b (c. a | c. b),mode=ramp,gradient=bwr,nbins=30,sat=.5, value=1.

      to color chains A and B with the Blue-White-Red gradient in 30 colors of equal 
      numbers of atoms in each color.
  """

  nbins=int(nbins)
  sat=float(sat)
  value=float(value)
# make sure sat and value are in the range 0-1.0
  sat = min(sat, 1.0)
  sat = max(sat, 0.0)
  value = min(value, 1.0)
  value = max(value, 0.0)

# make sure lowercase
  gradient.lower()
  mode.lower()

# Sanity checking
  if nbins == 1:
    print "\n     WARNING: You specified nbins=1, which doesn't make sense...resetting nbins=11\n"
    nbins=11

  if mode not in ('hist','ramp'):
    print "\n     WARNING: Unknown mode ",mode, "    ----->   Nothing done.\n"
    return
  elif gradient not in ('bgr','rgb','rainbow','reverserainbow','bwr','rwb',
                        'bmr','rmb','rw','wr','gw','wg','bw','wb','gy','yg','gray','grey','reversegray','reversegrey'):
    print "\n     WARNING: Unknown gradient: ",gradient, "    ----->   Nothing done.\n"
    return

  print "MODE, GRADIENT, NBINS:", mode,gradient, nbins

# get list of B-factors from selection
  m = cmd.get_model(selection)
  sel = []
  b_list = []

  if len(m.atom) == 0:
    print "Sorry, no atoms selected"

  else:
    if item == 'b':
      for i in range(len(m.atom)):
        b_list.append(m.atom[i].b)
    elif item == 'q':
      for i in range(len(m.atom)):
        b_list.append(m.atom[i].q)

    else:
      print "Not configured to work on item %s" % item
      return

    max_b = max(b_list)
    min_b = min(b_list)
    print "Minimum and Maximum B-values: ", min_b, max_b

    if mode == 'ramp':
      # color in bins of equal numbers of atoms
      b_list.sort()

      # subtract 0.1 from the lowest B in order to ensure that the single
      # atom with the lowest B value doesn't get omitted
#      b_list[0] = b_list[0] - 0.1

      bin_num = int(len(b_list)/nbins)
#      sel.append(selection + " and (b < " + str(b_list[bin_num]) + " or b = " + str(b_list[bin_num]) + ")")
      sel.append(selection + " and (%s < %4.4g" % (item,b_list[bin_num]) + " or %s = %4.4g" % (item,b_list[bin_num]) + ")")
      for j in range(1,nbins):
#        sel.append(selection + " and b > " + str(b_list[j*bin_num]))
        sel.append(selection + " and %s > %4.4g" % (item,b_list[j*bin_num]))
        #print "Color select: ",sel[j]

    elif mode == 'hist':

# check if minimum or maximum was specified and use the entered values
      if minimum != '':
        min_b = float(minimum)
      if maximum != '':
        max_b = float(maximum)
      # histogram:
      # color in bins of equal B-value ranges
      # subtract 0.1 from the lowest B in order to ensure that the single
      # atom with the lowest B value doesn't get omitted
      bin_width = (max_b - min_b)/nbins
      sel.append(selection + " and (%s < %4.4g" % (item,min_b + bin_width) + " or %s = %4.4g" % (item,min_b + bin_width) + ")")
      for j in range(1,nbins):
        sel.append(selection + " and %s > %4.4g" % (item,min_b + j*bin_width))
        #print "Color select: ",sel[j]

# call the function to create the gradient which returns a list of colours
    colours = make_gradient(sel,gradient,nbins,sat,value)

# do the colouring now
    for j in range(nbins):
      print "Color select: ",sel[j]
      cmd.color(colours[j],sel[j])

def color_q(selection="all",mode="hist",gradient="bgr",nbins=11,sat=1.,value=1.,minimum='',maximum=''):
  """
    
  USAGE

    color_q(selection,gradient,mode,nbins,sat,value) ='sel',
      gradient='bgr' or 'rgb' or 'bwr' or 'rwb' or 'bmr' or 'rmb' 
      'rw' or 'wr','gw' or 'wg' or 'bw' or 'wb' or 'gy' or 'yg' or 'gray' or 'reversegray'
      mode='hist' or 'ramp', q0=0.,q1=1.0,
      nbins=11, sat=1.0, value=1.0)

      This function allows coloring of a selection as a function of
      occupancy.  See color_b for details.
  """
  item='q'
  color_b(selection,item,mode,gradient,nbins,sat,value,minimum,maximum)

# function for creating the gradient
def make_gradient(sel,gradient,nbins,sat,value):
  if gradient == 'bgr' or gradient == 'rainbow':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # must append the str(sel[j]) to the color name so that it is unique 
      # for the selection
      coldesc.append('col' + str(j))
      # coldesc.append('col' + str(sel[j]) + str(j))

      # create colors using hsv scale (fractional) starting at blue(.6666667) 
      # through red(0.00000) in intervals of .6666667/(nbins -1) (the "nbins-1" 
      # ensures that the last color is, in fact, red (0)
      # rewrote this to use the colorsys module to convert hsv to rgb
      hsv = (colorsys.TWO_THIRD - colorsys.TWO_THIRD * float(j) / (nbins-1), sat, value)
      #convert to rgb and append to color list
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])
      # print "defined as ", str(sel[j])

  elif gradient == 'rgb' or gradient == 'reverserainbow':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # must append the str(sel[j]) to the color name so that it is unique 
      # for the selection
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))

      # create colors using hsv scale (fractional) starting at red(.00000) 
      # through blue(0.66667) in intervals of .6666667/(nbins -1) (the "nbins-1" 
      # ensures that the last color is, in fact, red (0)
      # rewrote this to use the colorsys module to convert hsv to rgb
      hsv = (colorsys.TWO_THIRD * float(j) / (nbins-1), sat, value)
      #convert to rgb and append to color list
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'bmr':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from blue through magenta to red
      rgb = [min(1.0, float(j)*2/(nbins-1)), 0.0, min(1.0, float(nbins-j-1)*2/(nbins-1))]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'rmb':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from red through magenta to blue
      rgb = [min(1.0, float(nbins-j-1)*2/(nbins-1)), 0.0, min(1.0, float(j)*2/(nbins-1))]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'rw':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from red through white
      rgb = [1.0, float(j)/(nbins-1), float(j)/(nbins-1)]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'wr':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from white through red 
      rgb = [1.0, float(nbins-j-1)/(nbins-1), float(nbins-j-1)/(nbins-1)]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'gw':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from green through white
      rgb = [float(j)/(nbins-1), 1.0, float(j)/(nbins-1)]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'wg':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from white through green 
      rgb = [float(nbins-j-1)/(nbins-1), 1.0, float(nbins-j-1)/(nbins-1)]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'bw':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from blue through white
      rgb = [float(j)/(nbins-1), float(j)/(nbins-1), 1.0 ]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'wr':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from white through blue 
      rgb = [float(nbins-j-1)/(nbins-1), float(nbins-j-1)/(nbins-1), 1.0]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'gy':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from green through yellow
      rgb = [float(j)/(nbins-1), 1.0, 0.]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'yg':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from green through yellow
      rgb = [float(nbins-j-1)/(nbins-1), 1.0, 0.]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'bwr':
    col=[]
    coldesc=[]
    for j in range(nbins/2):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from blue to white 
      rgb = [min(1.0, float(j)*2/(nbins-1)), min(1.0,float(j)*2/(nbins-1)), min(1.0, float(nbins-j-1)*2/(nbins-1))]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

    for j in range(nbins/2,nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from white to red
      rgb = [min(1.0, float(j)*2/(nbins-1)), min(1.0,float(nbins-j-1)*2/(nbins-1)), min(1.0, float(nbins-j-1)*2/(nbins-1))]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'rwb':
    col=[]
    coldesc=[]
    for j in range(nbins/2):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient from red to white
      rgb = [min(1.0, float(nbins-j-1)*2/(nbins-1)), min(1.0,float(j)*2/(nbins-1)), min(1.0, float(j)*2/(nbins-1))]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

    for j in range(nbins/2,nbins):
      coldesc.append('col' + str(j))
      # coldesc.append('col' + str(sel[j]) + str(j))
      # create colors in a gradient from white to blue
      rgb = [min(1.0, float(nbins-j-1)*2/(nbins-1)), min(1.0,float(nbins-j-1)*2/(nbins-1)), min(1.0, float(j)*2/(nbins-1))]

      # convert rgb to hsv,  modify saturation and value and convert back to rgb
      hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
      hsv[1] = hsv[1]*sat
      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'gray' or gradient == 'grey':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient of grays from "sat" to "value"

      hsv = [0, 0, sat + (value-sat)*float(j)/(nbins-1)]
#      hsv[1] = hsv[1]*sat
#      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])

  elif gradient == 'reversegray' or gradient == 'reversegrey':
    col=[]
    coldesc=[]
    for j in range(nbins):
      # coldesc.append('col' + str(sel[j]) + str(j))
      coldesc.append('col' + str(j))
      # create colors in a gradient of grays from "sat" to "value"

      hsv = [0, 0, value - (value-sat)*float(j)/(nbins-1)]
#      hsv[1] = hsv[1]*sat
#      hsv[2] = hsv[2]*value
      rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

      col.append(rgb)
      # cmd.set_color("col" + str(sel[j]) + str(j),col[j])
      cmd.set_color("col" + str(j),col[j])


# return the gradient as a list of colors named by their index (i.e. col0,col1,col2,col3,...)
  return coldesc

# allow calling without parentheses: color_hist_b [selection=], [mode= ],[gradient= ],[nbins= ]
cmd.extend("color_b",color_b)
cmd.extend("color_q",color_q)

