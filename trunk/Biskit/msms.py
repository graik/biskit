## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py
## DAG - substituted Numeric

## Class msms
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2012 Raik Gruenberg & Johan Leckner
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
##
## $Revision$
## last $Author$
## last $Date$
"""
Use MSMS to calculate surface info.

@note: This module is not any more the prefered module for
       surface area calculations. SurfaceRacer is the default
       module for this kind of calculations.
"""

import numpy as N
import Biskit.tools as T
import string
import os
import tempfile

## import Mslib

from Biskit import Executor, TemplateError, ExeConfigCache, ExeConfigError
## import Biskit.settings as S


class Pdb2xyzrnError( Exception ):
    pass

class Pdb2xyzrn( Executor ):
    """
    MSMS - pdb2xyzrn
    ================

    Run the awk script C{ pdb_to_xyzrn } to extract names and radii.
    This is nessesary to run msms. The radii and atoms depend on where
    in the sequence an atom resides. The radii used is determined
    by the information in the C{ atomtypenumbers } file provided
    with the MSMS application.

    Result
    ------
    ::
      r   - atom radii array
      n   - atom name list
    """

    def __init__( self, model, **kw ):
        """
        @param model: reference
        @type  model: PDBModel

        @param kw: additional key=value parameters for Executor:
        @type  kw: key=value pairs
        ::
          debug    - 0|1, keep all temporary files (default: 0)
          verbose  - 0|1, print progress messages to log (log != STDOUT)
          node     - str, host for calculation (None->local) NOT TESTED
                          (default: None)
          nice     - int, nice level (default: 0)
          log      - Biskit.LogFile, program log (None->STOUT) (default: None)

        """
        self.f_pdb = tempfile.mktemp('_pdb_to_xyzrn.pdb')

        ## gpdb_to_xyzrn have to be run i the local directory where
        ##   it resides. Otherwise it will not find the data file 
        ##   called "atmtypenumbers".
        if not os.path.exists( T.dataRoot() + '/msms/'):
            raise Pdb2xyzrnError, 'Cannot find msms directory. This should reside in ~biskit/Biskit/data/msms'

        fmsms =  T.dataRoot() + '/msms/'

        ## use biskit-version of pdb_to_xyzrn by default
        exe = ExeConfigCache.get('pdb2xyzrn', strict=1)
        try:
            exe.validate()
        except ExeConfigError:
            exe.bin = fmsms + '/gpdb_to_xyzrn'

        Executor.__init__( self, 'pdb2xyzrn', f_in=self.f_pdb,
                           cwd=fmsms, **kw )

        self.model = model


    def prepare( self ):
        """
        Overrides Executor method.
        """
        self.model.writePdb( self.f_pdb )


    def cleanup( self ):
        Executor.cleanup( self )

        if not self.debug:
            T.tryRemove( self.f_pdb )


    def parse_xyzrn( self, output ):
        """
        Extract xyz and r from output.

        @param output: STDOUT from pdb_to_xyzrn
        @type  output: [str]

        @return: r, n - array with atom radii, list with atom names
        @rtype: array, array
        """
        lines = output.split('\n')

        ## convert into lists that mslib will like
        xyzr = []
        n = []

        for l in lines:
            ## don't let the last empty line mess things up
            if len(l)!=0:
                l = string.split( l )
                xyzr += [ float(l[0]), float(l[1]), \
                          float(l[2]), float(l[3]) ]
                n += [ l[5] ]

        xyzr = N.reshape( xyzr, ( len(xyzr)/4, 4 ) )

        r = xyzr[:,3]

        ## xyz data is not passed on as it is not needed
        xyz = xyzr[:,:3]

        return r, n


    def isFailed( self ):
        """
        Overrides Executor method
        """
        return self.error is None 


    def finish( self ):
        """
        Overrides Executor method
        """
        Executor.finish( self )
        self.result = self.parse_xyzrn( self.output )


## class MSLIB_Error( Exception ):
##     pass


## class MSLIB:
##     """
##     MSLIB
##     =====
##     Calculate SAS (Solvent Accessible surface) and
##     SES (Solvent Excluded Surface) using MSLIB.

##     @note: This class has been retired. If you need to run MSMS
##            use the L[Biskit.MSMS} class insted. The current default class for
##            calculating MS and AS is otherwise L{Biskit.SurfaceRacer}

##     @bug: MSLIB fails when calling Mslib.MSMS()
##           The error seems to go back to Scientificpython:

##           ERROR::
##             File "/usr/tmp/python-8778q1e", line 160, in calc
##             surf = Mslib.MSMS( coords = xyz, radii = self.r )
##             File "Mslib/__init__.py", line 120, in __init__
##             msms.MOLSRF.__init__(self, name=name, coord=c, nat=nat, maxat=maxnat, names=atnames)
##             File "Mslib/msms.py", line 1832, in __init__
##             self.this = apply(msmsc.new_MOLSRF,_args,_kwargs)
##             ValueError: Failed to make a contiguous array of type 6
##     """

##     def __init__( self , model ): 

##         self.model = model

##         ## get radiia and name array
##         p2x = pdb2xyzrn(self.model)
##         self.r, self.n = p2x.run()


##     def calc( self, descr='no_name', probeRadius=1.5 ):
##         """
##         Run msms analytical surface calculation, i.e.::
##           msms -if xyzr.coord -af atoms.area -surface ases

##         @param descr: descriptive name
##         @type  descr: str
##         @param probeRadius: probe radius
##         @type  probeRadius: float

##         @return: out,  probe radii used, total SAS and SES
##                  sesList, array of lenght atoms with SES
##                         (Solvent Excluded Surface)
##                  sasList, array of lenght atoms with SAS
##                         (Solvent Accessible surface)
##                  atmList, list of atom names
##         @rtype: dict, array, array, [str]
##         """
##         xyz = self.model.xyz

##         print N.shape(xyz), N.shape(self.r)
##         print xyz[0], self.r[0]
##         print xyz[-1], self.r[-1]
##         ## run msms
##         surf = Mslib.MSMS( coords = xyz,
##                            radii = self.r, 
##                            atnames = self.n,
##                            name = descr,
##                            maxnat = len(xyz)+100 )

##         surf = Mslib.MSMS( coords = xyz, radii = self.r )

##         surf.compute( probe_radius = probeRadius, density=1.0 )
##         area = surf.compute_ses_area()

##         out = {}
##         if area == 0:
##             surf.write_ses_area( surfName )
##             out['ses'] = surf.sesr.fst.a_ses_area
##             out['sas'] = surf.sesr.fst.a_sas_area
##         else:
##             raise MSLIB_Error('compute_ses_area not successfull')

##         ## parse atom specific SAS (Solvent Accessible surface)
##         ## and SES (Solvent Excluded Surface)
##         surfName = surfName + '_0'
##         outList = open( surfName  , 'r').readlines()

##         sasList = []
##         sesList = []
##         atmList = []

##         for i in range( 1, len(outList) ):
##             line = string.split( outList[i] )
##             sesList += [ float( line[1] ) ]
##             sasList += [ float( line[2] ) ]
##             atmList += [ line[3] ]

##         ## polar and nonepolar surfaces
##         N_mask = N.transpose(atmList)[0] == 'N'
##         O_mask = N.transpose(atmList)[0] == 'O'
##         C_mask = N.transpose(atmList)[0] == 'C'
##         out['ses_polar'] = N.sum( sesList * O_mask ) + N.sum( sesList * N_mask )
##         out['ses_non-polar'] = N.sum( sesList * C_mask )
##         out['sas_polar'] = N.sum( sasList * O_mask ) + N.sum( sasList * N_mask )
##         out['sas_non-polar'] = N.sum( sasList * C_mask )

##         ## cleanup 
##         try:
##             os.remove( surfName )
##         except:
##             pass

##         return out, sesList, sasList, atmList


class MSMS_Error( Exception ):
    pass


class MSMS( Executor ):
    """
    MSMS
    ====

    Calculate SAS (Solvent Accessible surface) and
    SES (Solvent Excluded Surface) using the msms applicaton.

    Run msms analytical surface calculation i.e.::
       msms -if xyzr.coord -af atoms.area -surface ases

    Result
    ------
    ::
     out     - dictionary - probe radii used, total SAS and SES
     sesList - array of lenght atoms with SES (Solvent Excluded Surface)
     sasList - array of lenght atoms with SAS (Solvent Accessible surface)

    @note: The current default class for calculating MS and AS
           is L{Biskit.SurfaceRacer}.
    """


    def __init__( self, model, **kw ): 
        """
        @param model: PDBModel
        @type  model:

        @param kw: additional key=value parameters for Executor:
        @type  kw: key=value pairs
        ::
          debug    - 0|1, keep all temporary files (default: 0)
          verbose  - 0|1, print progress messages to log (log != STDOUT)
          node     - str, host for calculation (None->local) NOT TESTED
                          (default: None)
          nice     - int, nice level (default: 0)
          log      - Biskit.LogFile, program log (None->STOUT) (default: None)
        """
        self.f_xyzrn = tempfile.mktemp('_msms.xyzrn')

        ## output file from MSMS, will add .area exiension to file
        self.f_surf = tempfile.mktemp( )

        arg =' -surface ases -if %s  -af %s'%( self.f_xyzrn, self.f_surf )

        Executor.__init__( self, 'msms', args=arg, **kw )

        self.model = model.clone()


    def prepare( self ):
        """
        Write a xyzrn coordinate file to disc.
        Overrides Executor method.
        """
        ## get radiia and name array
        p2x = Pdb2xyzrn(self.model, verbose=self.verbose, debug=self.debug )
        r, n = p2x.run()

        xyz = self.model.xyz  
        xyzr = N.concatenate( ( xyz, N.transpose([r]) ) ,axis=1 )

        f = open( self.f_xyzrn, 'w' )
        i = 0
        for line in xyzr:
            f.write( str(line)[2:-1] + ' 1 ' + n[i] + '\n')
            i += 1
        f.close()


    def cleanup( self ):
        """
        Remove temp files.
        """
        Executor.cleanup( self )

        if not self.debug:
            T.tryRemove( self.f_xyzrn )
            T.tryRemove( self.f_surf + '.area' )


    def parse_msms( self ):
        """
        Parse the result file from a msms calculation.

        @return: probe radii used, total SAS and SES, 
                 array of lenght atoms with SES (Solvent Excluded Surface),
                 array of lenght atoms with SAS (Solvent Accessible surface),
                 and a list of atom names
        @rtype: dict, array, array, [str]
        """        
        ## get ASA, SES and probe radiii from stdout
        out = string.split( open( self.f_out ).readlines()[-3] )
        out = { 'probe':float(out[1]),
                'ses':float(out[-2]),
                'sas':float(out[-1]) }

        ## MSMS adds the extension .area to output file 
        ## over writing any existing extension
        surfName = self.f_surf + '.area'

        ## parse atom specific SAS (Solvent Accessible surface)
        ## and SES (Solvent Excluded Surface)
        outList = open( surfName  , 'r').readlines()

        sasList = []
        sesList = []
        atmList = []

        for i in range( 1, len(outList) ):
            line = string.split( outList[i] )

            if len(line)== 4:
                sesList += [ float( line[1] ) ]
                sasList += [ float( line[2] ) ]
                atmList += [ line[3] ]
            else:
                raise MSMS_Error, 'Wrong format in file %s'%surfName

        return out, sesList, sasList, atmList


    def isFailed( self ):
        """
        Overrides Executor method
        """
        return not self.error is None 


    def finish( self ):
        """
        Overrides Executor method
        """
        Executor.finish( self )
        self.result = self.parse_msms( )



#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test case"""

    TAGS = [ BT.OLD, BT.EXE ]  ## msms is superseeded by SurfaceRacer

    def test_msms(self):
        """msms test"""
        from Biskit import PDBModel

        if self.local: print 'Loading PDB...'
        self.m = PDBModel( T.testRoot() + '/lig/1A19.pdb' )
        self.m = self.m.compress( self.m.maskProtein() )

        if self.local: print 'Getting surfaces via the MSMS executable'
        self.ms = MSMS( self.m, debug=self.DEBUG, verbose=self.local )

        if self.local: print 'Running'
        out, sesList, sasList, atmList = self.ms.run()

        if self.local:
            print out
            print '\nResult from MSMS (first 5 lines): '
            print 'SES \tSAS \tAtom name'
            for i in range(5):
                print '%.2f \t%.2f \t%s'%(sesList[i], sasList[i], atmList[i])

            print 'MSMS done'

            globals().update( locals() )


        self.assertAlmostEqual( out['sas'] + out['ses'],
                                5085.1580000000004 + 4208.7389999999996, 8 )


if __name__ == '__main__':

    BT.localTest()

## ####################################
## ## Testing

## if __name__ == '__main__':

##     from Biskit import PDBModel

##     print "Loading PDB..."

##     m = PDBModel( T.testRoot() + '/lig/1A19.pdb' )
##     m = m.compress( m.maskProtein() )


## ##     #######
## ##     ## convert pdb file to coordinates, radii and names
## ##     print "Starting xpdb_to_xyzrn"
## ##     x = pdb2xyzrn( m, debug=0, verbose=1 )

## ##     print "Running"
## ##     r, n = x.run()


## ##     print '\nResult from pdb_to_xyzrn (first 5 lines): '
## ##     print 'radii \t\tnames'
## ##     for i in range(5):
## ##         print '%f \t%s'%(r[i], n[i])

##     #######
##     ## get surfaces via the MSMS executable

##     print "Starting MSMS"
##     a = MSMS( m, debug=1, verbose=1 )

##     print "Running"
##     out, sesList, sasList, atmList = a.run()

##     print out
##     print '\nResult from MSMS (first 5 lines): '
##     print 'SES \tSAS \tAtom name'
##     for i in range(5):
##         print '%.2f \t%.2f \t%s'%(sesList[i], sasList[i], atmList[i])

##     print 'MSMS done'


## ##     #######
## ##     ## get surfaces via MSLIB

## ##     a = MSLIB( m )
## ##     out, sesList, sasList, atmList = a.calc( m.xyz  )
## ##     print out
