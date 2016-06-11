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

import os
from Biskit.tools import *

time_offset = 40.0 ## (start time)

def get_rst_time( frst ):
    """extract time of last snapshot"""

    f = open( frst )
    f.readline()
    l = f.readline()

    f.close()

    t = float( l.split()[1] )

    return t

def rename_current_files( folder, current_time, exclude=[] ):
    """Rename *crd, *out, *vel, *rst to *_x_ps_.crd/out/vel/rst"""

    fs = os.listdir( folder )

    for f in fs:
        if f[-3:] in ['crd','out','vel','rst' ] and not f in exclude:

            ending = f[-3:]
            name = stripFilename( f )

            new_name = '%s/%s_%ips.%s' % \
                       (absfile( folder), name, current_time, ending)

            print "renaming %s -> %s" % (f, new_name)

            os.rename( folder + '/' + f, new_name )
    

def create_start_From_template( fTemplate, fout ):
    """get start.csh with right rst file"""

    pass

def adapt_inp( fold_inp, fnew_inp, current_time, old_total=0 ):
    """Parse inp, adapt step number to time remaining from total"""

    dt = 0

    sf = open( fold_inp ).read()
    sf = sf.replace('\n',',')
    sf = sf.replace(',,',',')
    sf = sf.replace('\n','')
    sf = sf.replace(' ', '')
    commands = sf.split(',')

    for c in commands:
        if c.find('=') != -1:
            param = c.split('=')[0]
            value = c.split('=')[1]

            if param == 'nstlim' and old_total==0:
                old_total = int( value )

            if param == 'dt':
                dt = float( value )
    
    if not( dt and old_total ):
        raise Exception('didnt find dt and/or nstlim option')


    new_total = int( ( old_total * dt - current_time ) / dt )

    fout = open( fnew_inp, 'w')
    
    for c in commands:

        if c.find('=') != -1 :
            fout.write('  ')

            if c.find('nstlim') != -1:
                fout.write('nstlim=%i,\n' % new_total)
            else:
                fout.write('%s,\n' % c )
        else:
            fout.write( '%s\n' % c )

    fout.close()
            

###MAIN###

o = cmdDict( {'f':'.', 't0':'40', 'rst':'sim.rst', 'inp':'sim.inp',
              'e':'eq.rst'} )

if len( sys.argv ) < 2:
    print \
    """
    Prepare the restart of a broken Amber MD run. Current *crd etc. are
    moved to oldName_TIMEps.* and the nstlim option in the input file
    is set to the number of steps remaining to the end of the MD.
    am_restartMD.py -f |folder| [ -t0 |time_offset| -tot |nstlim_total|
                       -rst |rst_file|
                       -inp |inp_file| -e |exclude_files_from_renaming| ]

       tot   - needed for 2nd restart, total number of MD steps (w/o restart)
       t0    - starting time in ps of this MD
       
    default:
    """
    for k,v in o.items():
        print "\t%s\t%s" % (k,v)

    sys.exit(0)

else:
    
    try:
        rst = absfile( o['f'] ) + '/' + o['rst']
    except:
        rst = absfile( o['rst'] )

    t0 = int( o['t0'] )
    tcurrent = get_rst_time( rst )

    finp = absfile( o['inp'] )

    exclude_files = toList( o['e'] )

    adapt_inp( finp, finp, tcurrent - t0, int( o.get('tot',0)) )

    rename_current_files( absfile( o['f'] ), tcurrent, exclude_files )
