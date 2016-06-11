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


"""
Calculate Contact Area Difference (CAD) with ICMBrowser.
"""

import tempfile, re

from Biskit import Executor, TemplateError
## import Biskit.settings as S
import Biskit.tools as T


class IcmCadError( Exception ):
    pass

class IcmCad( Executor ):
    """
    CAD calculation
    ===============
    Calculate ContactAreaDifference between a reference structure and one or
    many other structures. Described in Abagyan and Totrov, 1997.

    Example usage
    -------------
       >>> x = IcmCad( ref_model, [ m1, m2 ], verbose=1 )
       >>> result = x.run()

    IcmCad writes temporary pdb files, calls the icmbrowser program,
    and parses the result from the STDOUT pipe. The result values are then in
    x.result (if x is the IcmCad instance). As recommended, the value
    is multiplied by factor 1.8 which puts it into a range from 0 (identical)
    to 100 (completely different). This scaling is not recommended for
    docking solutions.

    @todo: Not speed-optimized!
             - it should be possible to pipe the PDB-files directly into
               icmbrowser instead of writing temporary files to disc
             - Abagyan and Totrov mention a c-library available -- the best
               solution would hence be to wrap this c-library in python...

    @note: Command configuration in: biskit/Biskit/data/defaults/exe_icmbrowser.dat
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
        @param ref_model: reference
        @type  ref_model: PDBModel
        @param models: structures to be compared with reference
        @type  models: [PDBModel]

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
        Executor.__init__( self, 'icmbrowser', template=self.inp_head, **kw )

        self.f_ref = tempfile.mktemp('_icmcad_ref.pdb')
        self.f_pdb = tempfile.mktemp('_icmcad_%i.pdb')

        self.ref_model = ref_model
        self.models = models

        if not isinstance( self.models, list ):
            self.models = [ self.models ]


    def __prepareModel( self, model, f_out ):
        """
        Prepare a model that ICM likes.
         - consecutive numbering of residues
         - no chain identifier

        @param model: model 
        @type  model: PDBModel
        @param f_out: name of pdb file to write
        @type  f_out: str        
        """
        model.renumberResidues()
        #for a in model.getAtoms():
            #a['chain_id']=''
        model['chain_id'] = [''] * len( model )
        model.writePdb( f_out )


    def prepare( self ):
        """
        Overrides Executor method.
        """
        self.__prepareModel( self.ref_model, self.f_ref )

        for i,m in enumerate( self.models ):

            self.__prepareModel( m, self.f_pdb % i )


    def cleanup( self ):
        """
        Tidy up the mess you created.
        """        
        Executor.cleanup( self )

        if not self.debug:
            T.tryRemove( self.f_ref )

            for i in range( len(self.models)):
                T.tryRemove( self.f_pdb % i )


    def parse_icm( self, output ):
        """
        Extract CAD value from icm output.

        @param output: STDOUT result of ICM run
        @type  output: [str]

        @return: ICM result
        @rtype: [str]

        @raise IcmCadError: if no result
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
        """
        Overrides Executor method
        """
        return self.error != '' or \
               not '3> 1.8*Cad' in self.output


    def fail( self ):
        """
        Overrides Executor method. Called when execution fails.
        """
        s = 'IcmBrowser failed. Please check the program output in the '+\
          'field `output` of this IcmCad instance (e.g. `print x.output`)!'
        self.log.add( s )

        raise IcmCadError, s

    def finish( self ):
        """
        Overrides Executor method
        """
        Executor.finish( self )
        self.result = self.parse_icm( self.output )


    def generateInp(self):
        """
        Has to be overridden to support head / segment / tail structure.
        Replace formatstr place holders in inp by fields of this class.

        @return: ICM input script
        @rtype: [str]

        @raise  TemplateError: if template not found
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


#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test class"""

    TAGS = [ BT.EXE ]

    def test_IcmCad( self ):
        """IcmCad test"""
        from Biskit import PDBModel

        if self.local: print 'Loading PDB...'

        self.f = T.testRoot() + '/lig/1A19.pdb'
        self.m1 = PDBModel(self.f)
        self.m1 = self.m1.compress( self.m1.maskProtein() )

        self.ms = []

        self.lig_traj = T.load( T.testRoot() + '/lig_pcr_00/traj.dat' )
        for m in self.lig_traj[:3]:
            m = m.compress( m.maskProtein() )
            self.ms.append(m)

        if self.local: print 'Starting ICM'
        self.x = IcmCad( self.m1, self.ms, debug=self.DEBUG,
                         verbose=self.local )

        if self.local:
            print 'Running'

        self.r = self.x.run()

        if self.local:
            print "Result: ", self.r

        self.assertEqual( self.r, [8.8603529999999999, 9.0315890000000003,
                                   8.5055429999999994] )



if __name__ == '__main__':

    BT.localTest()
