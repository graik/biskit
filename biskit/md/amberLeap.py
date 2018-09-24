##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2018 Raik Gruenberg & Johan Leckner
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


"""
Wrap the Amber tleap program.
"""

import tempfile
import os.path as P

from biskit.exe import Executor, TemplateError
import biskit.tools as T


class LeapError( Exception ):
    pass

class AmberLeap( Executor ):
    """
    Low-level wrapper for Amber tleap.

    AmberLeap expects a template for the leap input script which is then
    completed and passed to tleap. The manual command would like this::

      tleap -f leap.input

    The existence of a parm output file is used as success criterium.

    Usage::

        x = AmberLeap( 'leap_in.template' )
        x.run()

    However, typical leap scripts also load a PDB file and write
    top and crd files. Your template will thatswhy probably have placeholders
    for these. As usual, you can provide values for arbitrary place holders
    directly to the constructor::

        x = AmberLeap( 'leap_in.template',
                       leaprc='ff03',
                       out_parm='top.parm',
                       out_crd='0.crd',
                       in_pdb='1HPT_cleaned.pdb')
        x.run()


    Three special options are noteworthy:

      -leaprc -- the forcefield to be used,
                 will be completed to an existing leaprc file.
      -out_parm -- the name of the resulting topology file (default: top.parm)
      -out_crd -- the name of the resulting coordinate file (default: 0.crd)

    AmberParmBuilder uses AmberLeap and implements the gory details of PDB
    cleanup, S-S bond patching, chain break capping,
    and so on. It would be better to move this functionality directly
    into AmberLeap.

    .. seealso:: `biskit.Executor`, `biskit.AmberParmBuilder`
    """

    LEAPRC = 'leaprc.'             #: default leaprc file name
    LEAPRC_PATH = 'dat/leap/cmd/'  #: default leaprc location within $AMBERHOME

    def __init__( self,
                  template,
                  leaprc=None,
                  **kw ):
        """
        :param template: template for leap input file (file or string)
        :type  template: str
        :param leaprc: forcefield code (leaprc file ending, e.g. 'ff99')
                       OR leaprc file name (e.g, 'leaprc.f99')
                       OR leaprc path witin $AMBERHOME
                       OR leaprc path
                       default: take value from exe_tleap.dat
        :type  leaprc: str
        :param f_in:   complete leap input file -- existing or not, to be kept
                       (default: create from template and then discard)
        :type  f_in:   str
        :param f_out:  file for leap log output (default: discard)
        :type  f_out:  str

        :param out_parm: parm output file name (default: 'top.parm')
        :type  out_parm: str
        :param out_crd : coordinate output file name (default '0.crd' )
        :type  out_crd : str
        
        :param kw: additional key=value parameters for Executor:
        :type  kw: key=value pairs
        ::
          debug    - 0|1, keep all temporary files (default: 0)
          verbose  - 0|1, print progress messages to log (log != STDOUT)
          node     - str, host for calculation (None->local) NOT TESTED
                          (default: None)
          nice     - int, nice level (default: 0)
          log      - biskit.LogFile, program log (None->STOUT) (default: None)
        """
        # override some Executor defaults unless they are freshly given
        kw['catchout'] = kw.get('catchout',0)

        Executor.__init__( self, 'tleap', template=template, **kw )

        self.args = '-f %s' % self.f_in

        self.leaprc = self.findLeaprc( leaprc or self.exe.leaprc )

        ## set some defaults that may or may not have been specified
        self.out_parm = kw.get('out_parm', 'top.parm' )
        self.out_crd  = kw.get('out_crd',  '0.crd')
        

    def findLeaprc( self, ff ):
        """
        Guess full path to an existing leaprc file name or fragment of it.
        We first take the file as is, then look in AMBERHOME, then look
        in AMBERHOME/dat/leap/cmd.
        :param ff: forcefield code (leaprc file ending, e.g. 'ff99', 'ff01')
                   OR leaprc file name (e.g, 'leaprc.f99')
                   OR leaprc path witin $AMBERHOME
                   OR leaprc path
        :type  ff: str
        :return: full path to an existing leaprc file
        :rtype: str
        :raise: LeapError, if no existing file can be found
        """
        if P.exists( T.absfile(ff) ):
            return T.absfile(ff)

        amberhome = self.exe.env['AMBERHOME']

        ## interpret ff as file ending ('ff99', 'ff03', etc)
        r = P.join( amberhome, self.LEAPRC_PATH, self.LEAPRC + ff )
        if P.exists( r ):
            return r

        ## interpret ff as file ending ('ff99', 'ff03', etc pointing to 'old' forcefield)
        r = P.join( amberhome, self.LEAPRC_PATH, 'oldff', self.LEAPRC + ff )
        if P.exists( r ):
            return r

        ## interpret ff as file name (e.g. 'leaprc.ff99')
        r = P.join( amberhome, self.LEAPRC_PATH, ff )
        if P.exists( r ):
            return r

        ## interpret ff as file name (e.g. 'leaprc.ff99') for 'old' forcefield
        r = P.join( amberhome, self.LEAPRC_PATH, 'oldff', ff )
        if P.exists( r ):
            return r

        ## interpret ff as path within AMBERHOME ('dat/leap/cmd/leaprc.f99')
        r = P.join( amberhome, ff )
        if P.exists( r ):
            return r

        raise  LeapError('Could not find Forcefield definition %s. ' % ff +\
                  'Check exe_tleap.dat or provide an explicit leaprc parameter!')
            

    def isfailed( self ):
        return not os.path.exists( self.out_parm )
    
    def cleanup(self):
        if not self.debug:
            T.tryRemove( P.join(self.cwd or '', 'leap.log'))
        
        super().cleanup()


#############
## TESTING ##
import biskit.test as BT
import tempfile
import biskit.tools as T

class Test( BT.BiskitTest ):
    """Test AmberParmBuilder"""

    TAGS = [ BT.EXE ]

    def prepare(self):
        root = T.testRoot() + '/amber/'
        self.fpdb = T.testRoot() + '/amber/1HPT_0.pdb'
        self.fparm = tempfile.mktemp('.parm', 'top_')
        self.flog  = tempfile.mktemp('.log', 'leap_')
        self.template = T.dataRoot() + '/amber/leap/solvate_box.leap'
        
    def cleanUp(self):
        import os.path as osp
        T.tryRemove( self.fparm )
        T.tryRemove( self.flog )
        

    def test_AmberLeap_findLeaprc(self):
        """AmberLeap.findLeaprc test"""
        self.x = AmberLeap(self.template)
        target = P.join(self.x.exe.env['AMBERHOME'], 'dat/leap/cmd/oldff/leaprc.ff10')

        if self.local:
            self.log.add( '\ntarget leaprc: %s' % target )
        
        self.assertEqual( self.x.findLeaprc('ff10'), target )
        self.assertEqual( self.x.findLeaprc('leaprc.ff10'), target )
        self.assertEqual( self.x.findLeaprc('dat/leap/cmd/oldff/leaprc.ff10'),target)
        self.assertEqual( self.x.findLeaprc( target ), target )

    def test_AmberLeap_run( self ):
        """AmberLeap.run test"""
        self.x = AmberLeap( self.template,
                            leaprc='protein.ff03.r1',
                            leap_out=self.flog,
                            fmod='', fprep='', ss_bonds='', box=12.5,
                            in_pdb=self.fpdb,
                            debug=self.DEBUG)
        self.run()
        
        

if __name__ == '__main__':

    BT.localTest()
