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
"""Calculate Contact Area Difference (CAD) with ICMBrowser."""

import tempfile, re

from Biskit import Executor, TemplateError
import Biskit.settings as S


class IcmCadError( Exception ):
    pass

class IcmCad( Executor ):
    """
    Calculate ContactAreaDifference between a reference structure and one or
    many other structures.
    Described in Abagyan and Totrov, 1997.

    Example usage:
    
    x = IcmCad( ref_model, [ m1, m2 ], verbose=1 )
    result = x.run()

    IcmCad writes temporary pdb files, calls the icmbrowser program,
    and parses the result from the STDOUT pipe. The result values are then in
    x.result (if x is the IcmCad instance). As recommended, the value
    is multiplied by factor 1.8 which puts it into a range from 0 (identical)
    to 100 (completely different). This scaling is not recommended for
    docking solutions.

    ToDo: Not speed-optimized!
    a) it should be possible to pipe the PDB-files directly into icmbrowser
       instead of writing temporary files to disc
    b) Abagyan and Totrov mention a c-library available -- the best solution
       would hence be to wrap this c-library in python...

    Command configuration in: biskit/external/defaults/exe_icmbrowser.dat
    """

    inp_head = \
"""
call _startup
read pdb unix cat \"%(f_ref)s\"
copy a_ \"mol1\"
"""
    inp_tail = "quit\n"

    inp_calc = \
"""
read pdb unix cat \"%(f_pdb)s\"
copy a_ \"mol2\" delete
1.8*Cad( a_mol1./* a_mol2./* )
"""

    def __init__( self, ref_model, models, **kw ):
        """
        ref_model - PDBModel, reference
        models    - [ PDBModel ], structures to be compared with reference

        ... and additional key=value parameters for Executor:
        debug    - 0|1, keep all temporary files                       [0]
        verbose  - 0|1, print progress messages to log     [log != STDOUT]
        node     - str, host for calculation (None->local) NOT TESTED [None]
        nice     - int, nice level                                     [0]
        log      - Biskit.LogFile, program log (None->STOUT)        [None]
        """
        Executor.__init__( self, 'icmbrowser', template=self.inp_head, **kw )

        self.f_ref = tempfile.mktemp('_icmcad_ref.pdb')
        self.f_pdb = tempfile.mktemp('_icmcad_%i.pdb')

        self.ref_model = ref_model
        self.models = models

        if not isinstance( self.models, list ):
            self.models = [ self.models ]


    def __prepareModel( self, model, f_out ):
        model.renumberResidues()
        for a in model.getAtoms():
            a['chain_id']=''
        model.writePdb( f_out )

    def prepare( self ):
        """Overrides Executor method."""
        self.__prepareModel( self.ref_model, self.f_ref )

        for i,m in enumerate( self.models ):
            
            self.__prepareModel( m, self.f_pdb % i )


    def cleanup( self ):
        Executor.cleanup( self )

        if not self.debug:
            T.tryRemove( self.f_ref )

            for i in range( len(self.models)):
                T.tryRemove( self.f_pdb % i )


    def parse_icm( self, output ):
        """
        Extract CAD value from icm output.
        lines - [str], STDOUT result of ICM run
        """
        lines = output.split('\n')
        if  len(lines) == 0:
            raise IcmCadError, 'no ICM result'

        r = []

        for i, l in enumerate(lines):  
            if re.match('.+1.8\*Cad.+', l ):
                r.append( float( lines[i+1].strip() ) )

        return r

    def isFailed( self ):
        """Overrides Executor method"""
        return not self.error is None 

    def finish( self ):
        """Overrides Executor method"""
        Executor.finish( self )
        self.result = self.parse_icm( self.output )

    
    def generateInp(self):
        """
        Has to be overridden to support head / segment / tail structure.
        Replace formatstr place holders in inp by fields of this class.
        inp  - ignored
        fout - ignored
        -> ICM input script 
        !! raise TemplateError
        """
        try:
            s = self.inp_head % self.__dict__

            for i,m in enumerate( self.models ):
                s += self.inp_calc % { 'f_pdb' : self.f_pdb % i }

            s += self.inp_tail

            return s
        
        except KeyError, why:
            s =  "Unknown option/place holder in template file."
            s += "\n  Template asked for a option called " + str( why[0] )
            raise TemplateError, s
        
        except Exception, why:
            s =  "Error while creating template file."
            s += "\n  why: " + str( why )
            s += "\n  Error:\n  " + T.lastError()
            raise TemplateError, s



if __name__ == '__main__':

    from Biskit import PDBModel
    import Biskit.tools as T
    import glob

    print "Loading PDBs..."

    m1 = PDBModel( T.testRoot()+'/lig_pc2_00/pdb/1A19_1_0.pdb.gz' )
    m1 = m1.compress( m1.maskProtein() )

    ms = []

    for f in glob.glob( T.testRoot()+'/lig_pc2_00/pdb/*_1_*pdb.gz' )[:3]:
        m = PDBModel(f)
        m = m.compress( m.maskProtein() )
        ms.append(m)
        T.flushPrint('*')

    print "Starting ICM"

    x = IcmCad( m1, ms, debug=0, verbose=1 )

    print "Running"
    r = x.run()

    print "Result: ", r
    
