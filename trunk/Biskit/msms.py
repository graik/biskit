## Class msms
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
##
## $Revision$
## last $Author$
## last $Date$

import Numeric as N
import string
import os
import tempfile

import Mslib

from Biskit import Executor, TemplateError
import Biskit.settings as S


class pdb2xyzrnError( Exception ):
    pass

class pdb2xyzrn( Executor ):
    """
    Run the awk script pdb_to_xyzrn to extract names and radii.
       nessesary to run msms. The radii od atoms depend on where
       in the sequence an atom resides. The radii used is determined
       by the information in the 'atomtypenumbers' file provided
       with together with the MSMS application.
       
    ->  r   - atom radii array
        n   - atom name list
    """
    
    def __init__( self, model, **kw ):
        """
        model - PDBModel, reference

        ... and additional key=value parameters for Executor:
        debug    - 0|1, keep all temporary files                       [0]
        verbose  - 0|1, print progress messages to log     [log != STDOUT]
        node     - str, host for calculation (None->local) NOT TESTED [None]
        nice     - int, nice level                                     [0]
        log      - Biskit.LogFile, program log (None->STOUT)        [None]
        """
        self.f_pdb = tempfile.mktemp('_pdb_to_xyzrn.pdb')
        
        Executor.__init__( self, 'pdb2xyzrn', f_in=self.f_pdb, **kw )

        self.model = model


    def prepare( self ):
        """Overrides Executor method."""
        self.model.writePdb( self.f_pdb )


    def cleanup( self ):
        Executor.cleanup( self )

        if not self.debug:
            T.tryRemove( self.f_pdb )

            
    def parse_xyzrn( self, output ):
        """
        Extract xyz and r from output.
        output - [str], STDOUT from pdb_to_xyzrn

        -> r - array with atom radii
           n - list with atom names
        """
        lines = output.split('\n')
        ## convert into lists that mslib will like
        xyzr = []
        r = []
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
        """Overrides Executor method"""
        return not self.error is None 


    def finish( self ):
        """Overrides Executor method"""
        Executor.finish( self )
        self.result = self.parse_xyzrn( self.output )

        
class MSLIB_Error( Exception ):
    pass


class MSLIB:
    """
    Calculate SAS (Solvent Accessible surface) and
    SES (Solvent Excluded Surface) using msms.

    NOTE: MSLIB fails when calling Mslib.MSMS()
          The error seems to go back to Scientificpython.
          ERROR:
            File "/usr/tmp/python-8778q1e", line 160, in calc
            surf = Mslib.MSMS( coords = xyz, radii = self.r )
            File "Mslib/__init__.py", line 120, in __init__
            msms.MOLSRF.__init__(self, name=name, coord=c, nat=nat, maxat=maxnat, names=atnames)
            File "Mslib/msms.py", line 1832, in __init__
            self.this = apply(msmsc.new_MOLSRF,_args,_kwargs)
            ValueError: Failed to make a contiguous array of type 6
    """

    def __init__( self , model ): 
     
        self.model = model
 
        ## get radiia and name array
        p2x = pdb2xyzrn(self.model)
        self.r, self.n = p2x.run()
      

    def calc( self, descr='no_name', probeRadius=1.5 ):
        """
        Run msms analytical surface calculation
        i.e. msms -if xyzr.coord -af atoms.area -surface ases

        xyzrName - path and name of xyzr coordinate file

        ->  out     - dictionary - probe radii used, total SAS and SES
            sesList - array of lenght atoms with SES
                        (Solvent Excluded Surface)
            sasList - array of lenght atoms with SAS
                        (Solvent Accessible surface)
        """
        xyz = self.model.xyz

        print N.shape(xyz), N.shape(self.r)
        print xyz[0], self.r[0]
        print xyz[-1], self.r[-1]
        ## run msms
        surf = Mslib.MSMS( coords = xyz,
                           radii = self.r, 
                           atnames = self.n,
                           name = descr,
                           maxnat = len(xyz)+100 )

        surf = Mslib.MSMS( coords = xyz, radii = self.r )

        surf.compute( probe_radius = probeRadius, density=1.0 )
        area = surf.compute_ses_area()

        out = {}
        if area == 0:
            surf.write_ses_area( surfName )
            out['ses'] = surf.sesr.fst.a_ses_area
            out['sas'] = surf.sesr.fst.a_sas_area
        else:
            raise MSLIB_Error('compute_ses_area not successfull')
        
        ## parse atom specific SAS (Solvent Accessible surface)
        ## and SES (Solvent Excluded Surface)
        surfName = surfName + '_0'
        outList = open( surfName  , 'r').readlines()

        sasList = []
        sesList = []
        atmList = []

        for i in range( 1, len(outList) ):
            line = string.split( outList[i] )
            sesList += [ float( line[1] ) ]
            sasList += [ float( line[2] ) ]
            atmList += [ line[3] ]

        ## polar and nonepolar surfaces
        N_mask = N.transpose(atmList)[0] == 'N'
        O_mask = N.transpose(atmList)[0] == 'O'
        C_mask = N.transpose(atmList)[0] == 'C'
        out['ses_polar'] = N.sum( sesList * O_mask ) + N.sum( sesList * N_mask )
        out['ses_non-polar'] = N.sum( sesList * C_mask )
        out['sas_polar'] = N.sum( sasList * O_mask ) + N.sum( sasList * N_mask )
        out['sas_non-polar'] = N.sum( sasList * C_mask )
                
        ## cleanup 
        try:
            os.remove( surfName )
        except:
            pass

        return out, sesList, sasList, atmList

 
class MSMS_Error( Exception ):
    pass


class MSMS( Executor ):
    """
    Calculate SAS (Solvent Accessible surface) and
    SES (Solvent Excluded Surface) using msms.

    Run msms analytical surface calculation
    i.e. msms -if xyzr.coord -af atoms.area -surface ases

    ->  out     - dictionary - probe radii used, total SAS and SES
        sesList - array of lenght atoms with SES (Solvent Excluded Surface)
        sasList - array of lenght atoms with SAS (Solvent Accessible surface)
    """
                  

    def __init__( self, model, **kw ): 
        """
        model - PDBModel

        ... and additional key=value parameters for Executor:
        debug    - 0|1, keep all temporary files                       [0]
        verbose  - 0|1, print progress messages to log     [log != STDOUT]
        node     - str, host for calculation (None->local) NOT TESTED [None]
        nice     - int, nice level                                     [0]
        log      - Biskit.LogFile, program log (None->STOUT)        [None]
        """
        self.f_xyzrn = tempfile.mktemp('_msms.xyzrn')

        ## output file from MSMS, will add .area exiension to file
        self.f_surf = tempfile.mktemp( )
        
        arg =' -surface ases -if %s  -af %s'%( self.f_xyzrn, self.f_surf )
        
        Executor.__init__( self, 'msms', args=arg, **kw )

        self.model = model.clone( deepcopy=1 )

            
    def prepare( self ):
        """
        Overrides Executor method.
        Write a xyzrn coordinate file to disk
        """
        ## get radiia and name array
        p2x = pdb2xyzrn(self.model)
        r, n = p2x.run()
        
        xyz = self.model.xyz  
        xyzr = N.concatenate(  ( xyz, N.transpose([r]) ) ,axis=1 )

        f = open( self.f_xyzrn, 'w' )
        i = 0
        for line in xyzr:
            f.write( str(line)[2:-1] + ' 1 ' + n[i] + '\n')
            i += 1
        f.close()


    def cleanup( self ):
        Executor.cleanup( self )

        if not self.debug:
            T.tryRemove( self.f_xyzrn )
            T.tryRemove( self.f_surf + '.area' )

                
    def parse_msms( self ):
        """
        ->  out     - dictionary - probe radii used, total SAS and SES
            sesList - array of lenght atoms with SES
                       (Solvent Excluded Surface)
            sasList - array of lenght atoms with SAS
                       (Solvent Accessible surface)
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
        """Overrides Executor method"""
        return not self.error is None 


    def finish( self ):
        """Overrides Executor method"""
        Executor.finish( self )
        self.result = self.parse_msms( )
        

####################################
## Testing

if __name__ == '__main__':

    from Biskit import PDBModel
    import Biskit.tools as T
    import glob
    import mathUtils as MA

    print "Loading PDB..."

    f = glob.glob( T.testRoot()+'/lig_pc2_00/pdb/*_1_*pdb.gz' )[1]
    m = PDBModel(f)
    m = m.compress( m.maskProtein() )

##     #######
##     ## convert pdb file to coordinates, radii and names
##     print "Starting xpdb_to_xyzrn"
##     x = pdb2xyzrn( m, debug=0, verbose=1 )

##     print "Running"
##     r, n = x.run()


##     print '\nResult from pdb_to_xyzrn (first 5 lines): '
##     print 'radii \t\tnames'
##     for i in range(5):
##         print '%f \t%s'%(r[i], n[i])

    #######
    ## get surfaces via the MSMS executable
    
    print "Starting MSMS"
    a = MSMS( m, debug=1, verbose=1 )

    print "Running"
    out, sesList, sasList, atmList = a.run()
    
    print out
    print '\nResult from MSMS (first 5 lines): '
    print 'SES \tSAS \tAtom name'
    for i in range(5):
        print '%.2f \t%.2f \t%s'%(sesList[i], sasList[i], atmList[i])

    print 'MSMS done'


##     #######
##     ## get surfaces via MSLIB

##     a = MSLIB( m )
##     out, sesList, sasList, atmList = a.calc( m.xyz  )
##     print out
    
