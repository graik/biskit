##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2016 Raik Gruenberg & Johan Leckner

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

"""work in progress!"""

from Biskit import Xplorer, XplorerError
from Biskit import PDBModel, PCRModel
import Biskit.tools as t
from Biskit.Dock.ComplexEvolving import ComplexEvolving

import tempfile
import os
import copy
import shutil

class XRefineError( XplorerError ):
    pass

class XRefineComplex( Xplorer ):
    """
    Run XPlor refinement of a single docking solution.
    """

    def __init__( self, com, inp=None, **params ):
        """
        com  - Complex, to be refined
        inp  - alternative input template (also override path='..')

        passed on to Xplorer:
        xout    - str, file name of xplor log file (None->discard) [None]
        bin     - str, file name of xplor binary                  [xplor]
        node    - str, host for calculation (None->local)          [None]
        nice    - int, nice level                                     [0]
        log     - Biskit.LogFile, program log (None->STOUT)        [None]
        debug   - 0|1, keep all temporary files                       [0]
        verbose - 0|1, print xplor command to log (0-> log != STDOUT) [0]

        [ optional key=value pairs for inp file ], among them:
        path    - str, location of Xplor helper scripts [inpRefine/ in biskit]
        com_psf - str, output file name for complex psf
        com_out - str, output file name for refined complex pdb
        """
        self.com = com
        
        self.rec_psf = com.rec().getPsfFile()
        self.lig_psf = com.lig().getPsfFile()

        recCode = com.rec().getPdbCode()
        ligCode = com.lig().getPdbCode()

        ## tempfile.tempdir = settings.tempDirShared ## DEBUGGING

        self.rec_in = tempfile.mktemp( recCode + ".pdb" )
        self.lig_in = tempfile.mktemp( ligCode + ".pdb" )

        self.lig_out = tempfile.mktemp( "lig_out.pdb" )
        self.rec_out = tempfile.mktemp( "rec_out.pdb" )

        self.com_psf = params.get('com_psf', tempfile.mktemp( "com_out.psf" ) )
        self.com_out = params.get('com_out', tempfile.mktemp( "com_out.pdb" ) )

        self.ene_out = tempfile.mktemp("energies.out")

        ## location of additional XPlor scripts, required by refine.inp
        self.path = t.dataRoot() + '/xplor/inpRefine/'

        if len( self.path ) > 80-16:
            self.__localizeParams( params )

        assert len( self.path ) <= 80-16, '%s too long, ' % self.path +\
               'Xplor has a 80 character limit for pathnames.'

        self.inp_template = inp or self.path + 'refine.inp'

        self.result = None

        Xplorer.__init__(self, self.inp_template, **params )


    def version( self ):
        return Xplorer.version(self) + '; XRefineComplex $Revision$'


    def __localizeParams( self, params ):
        """
        Workaround for xplor limit of 80 characters. Copy the parameter
        folder into the working folder so that pathnames in the input file
        can stay below 80 characters.
        """
        params['cwd'] = params.get('cwd', tempfile.tempdir )
        self.path = os.path.join(params['cwd'], 'inpRefine') + '/'

        try:
            t.tryRemove( self.path, tree=True )
            shutil.copytree( t.dataRoot() + '/xplor/inpRefine',
                             self.path )
        except OSError, why:
            raise XRefineError('cannot copy parameter folder to working'+\
                               ' directory %s. ' % params['cwd'] + str(why) )


    def prepare( self ):
        
        self.com.rec().writePdb( self.rec_in )
        self.com.lig().writePdb( self.lig_in )

        
    def cleanup( self ):

        Xplorer.cleanup( self )

        if not self.debug:

            t.tryRemove( self.rec_in )
            t.tryRemove( self.lig_in )

            t.tryRemove( self.rec_out )
            t.tryRemove( self.lig_out )

            t.tryRemove( self.com_psf )
            t.tryRemove( self.com_out )

            t.tryRemove( self.ene_out )

            if self.cwd in self.path:
                t.tryRemove( self.path, tree=True )


    def readEnergies( self ):
        """
        Extract energies from xplor output file
        -> { 'str' : float }
        !! XRefineError
        """
        try:
            lines = open( self.ene_out, 'r' ).readlines()
            keys = lines[0].split()

            r = {}
            for k in keys:
                r[ k ] = []

            for l in lines[1:]:
                values = [ float( v ) for v in l.split() ]

                for k,v in zip( keys, values ):
                    r[ k ] += [ v ]

            return r
            
        except IOError, why:
            raise XRefineError, 'Cannot open xplor energy file %s' \
                  % self.ene_out
        except Exception, why:
            raise XRefineError, 'Cannot read energies from file %s' \
                  % self.ene_out + '\nReason: ' + t.lastError()


    def cloneAtoms( self, template_model, model ):
        """
        Make model look exactly like template_model but don't touch the
        coordinates.
        template_model - PDBModel, use these residue names, chain_id, segid..
        model          - PDBModel, same atom content as template in same order
        -> PDBModel, coordinates from model & atom dictionaries from template
        """

        i_model, i_templ = model.compareAtoms( template_model )
        
        if i_templ != template_model.atomRange():
            raise XRefineError,\
                  'cannot equalize two models with different atom content.'

        model = model.take( i_model )
        
        m = template_model.clone()
        m.setXyz( model.getXyz() )

        assert m.atomNames() == template_model.atomNames()

        return m


    def finish( self ):

        rec = self.cloneAtoms( self.com.rec(), PDBModel(source=self.rec_out) )

        lig = self.cloneAtoms( self.com.lig(), PDBModel(source=self.lig_out) )

        if len( rec ) and len( lig ):
            self.result = ComplexEvolving( rec, lig, self.com,
                                           info={'comment':'refined 1'} )

            self.result['refine_energies'] = self.readEnergies()
        else:
            self.result = None

        
    def fail( self ):
        Xplorer.fail( self )

#### TEST #######

if __name__ == '__main__':

##     cl = t.load( '~/interfaces/c11/dock_xray/hex01/complexes_cont.cl' )
    cl = t.load( t.testRoot() + '/dock/hex/complexes.cl')

    c = cl[0]

    x = XRefineComplex(c, ##"~/data/tmp",
                       inp=\
                       t.dataRoot() + '/xplor/inpRefine/refine.inp',
##                        xout = 'test_xplor.log',
                       debug=1,
                       verbose=1,
                       pipes=0)

    x.run()
