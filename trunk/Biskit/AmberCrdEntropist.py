##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner
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
## last $Date$
## $Revision$

import os.path as osp
import re
import tempfile

import Biskit.tools as T
import Biskit.settings as settings
from Biskit.Errors import BiskitError
from Biskit.Executor import Executor

class EntropistError( BiskitError ):
    pass

class AmberCrdEntropist( Executor ):
    """
    Run ptraj entropy analysis on existing amber trajectory (.crd) and
    topology (.parm7) file.
    """

    ## template for ptraj input script
    ptraj_script = """
    trajin %(f_crd)s %(start)i %(stop)i %(step)i
    matrix mwcovar name mwc
    analyze matrix mwc vecs 0 thermo
    """

    def __init__( self, f_parm, f_crd, f_template=None,
                  s=0, e=None, step=1, **kw ):
        """
        f_parm     - str, amber topology file
        f_crd      - str, amber trajectory file
        f_template - str, alternative ptraj input template [default]
        s,e        - int, start and stop frame [0, None .. all ]
        step       - int, frame offset         [1 .. no offset ]

        ... and key=value parameters for Executor:
        f_out    - str, target name for ptraj output file  [None..discard]
        node     - str, host for calculation (None->local)          [None]
        nice     - int, nice level                                     [0]
        log      - Biskit.LogFile, program log (None->STOUT)        [None]
        debug    - 0|1, keep all temporary files                       [0]
        verbose  - 0|1, print progress messages to log     [log != STDOUT]
        """
        template = f_template or self.ptraj_script

        Executor.__init__( self, 'ptraj', template=template, push_inp=0, **kw )

        self.f_parm = T.absfile( f_parm )
        self.f_crd  = T.absfile( f_crd )

        self.start = s
        self.stop  = e or 10000000
        self.step  = step

        ## default result
        self.result = { 'T':None, 'mass':None, 'vibes':None,
                        'S_total':None, 'S_trans':None, 'S_rot':None,
                        'S_vibes':None,
                        'contributions':None, 'nframes':None,
                        'version':self.version(), 'node':self.node }

                           
    def version( self ):
        return 'AmberCrdEntropist $Revision$'

    def command( self ):
        return "%s %s %s" % (self.exe.bin, self.f_parm, self.f_in)


    def __tryMatch( self, regex, str, integer=0 ):
        m = regex.search( str )
        if m:
            return round( float( m.groups()[0] ), 3 )
        return None


    def finish( self ):
        self.result = self.parsePtrajResult( self.f_out )


    def isFailed( self ):
        """Detect whether external program failed, override!"""
        if not osp.exists( self.f_out ):
            return 1

        return 0

    def parsePtrajResult( self, f_out ):

        ## regular expressions for parsing of ptraj output
        re_thermo= re.compile('- Thermochemistry -')
        re_T     = re.compile('^\s*temperature\s*(\d+\.\d+)\s*kelvin')
        re_mass  = re.compile('^\s*molecular mass\D*(\d+\.\d+)\s*amu')
        re_vibes = re.compile('^\s*Warning--\s*(\d+)\s*vibrations have low')
        re_table = re.compile('^-{80}$')
        re_value = re.compile('(-*\d+\.\d+)$')
        re_nsets = re.compile('Successfully read in (\d+) sets')

        r = self.result

        try:
            lines = open( f_out, 'r' ).readlines()
            lines.reverse()
            l = lines.pop()

            while not re_thermo.search(l):
                l = lines.pop()

            while not re_table.search(l):
                r['T']    = r['T']     or self.__tryMatch( re_T, l )
                r['mass'] = r['mass']  or self.__tryMatch( re_mass, l )
                r['vibes']= r['vibes'] or self.__tryMatch( re_vibes,l )
                l = lines.pop()

            l = lines.pop()
            v = []
            while l != '\n':
                v += [ self.__tryMatch( re_value, l ) ]
                l = lines.pop()

            r['S_total'], r['S_trans'], r['S_rot'], r['S_vibes'] = tuple(v[:4])
            r['contributions'] = v[4:]

            while len( lines ) > 0:
                r['nframes'] = r['nframes'] or self.__tryMatch( re_nsets, l )
                l = lines.pop()

        except IndexError:
            raise EntropistError, 'unexpected end of ptraj output file.'

        return r

if __name__ == '__main__':

    print "Setting up"

    f = T.testRoot() + '/Amber/AmberCrdEntropist/' 
    
    e = AmberCrdEntropist( f + 'lig_traj.parm',
                           f + 'lig_traj.crd', debug=0, verbose=1)

    print "Running"

    e.run()
