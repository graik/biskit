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

## last $Author$
## last $Date$
## $Revision$

"""
Triangulate a protein-protein interface with Intervor.
"""

import tempfile, re
import numpy as N
import os.path as osp
import glob

import Biskit as B
import Biskit.tools as T
from Biskit import Executor, TemplateError


class IntervorError( Exception ):
    pass

class Intervor( Executor ):
    """
    Triangulate a protein-protein interface with Intervor.
    The input protein must contain waters.
    
    Use::

      x = Intervor( model, [0,1] )
      x.run()

    ... where the [0,1] designates the receptor chains, the remaining
    protein chains are considered ligand, the remaining chains water.

    The parser creates several atom profiles in model:

        * facets     ... raw info of all facets the atom is involved in

        * n_facets   ... number of facets the atom is involved in
        * facet_area ... total area of this atom's facets
        * so_min     ... lowest shelling order
        * so_max     ... highest shelling order 

    For quick visualization use::

        pymoler = x.visualize( profile='so_max', xwat=[3], wat=[4] )
        pymoler.show()  ## needs pymol installed

    Which gives you a pymol instance with rec, lig, crystal waters
    (chain 3), and remaining waters (chain 4) as separate groups and
    color-coded by 'so_max'.

    """
    
    #: reg expressions for parsing
    RE_SECTION_SO = re.compile( '^Facet count' )
    
    RE_facet = re.compile( '^(\d+) *\| *(\d+) *([\w\d]+) *(\w+) *(\w) *(\d+) *\| *(\d+) *([\w\d]+) *(\w+) *(\w) *(\d+) *\| *(\d+) +(\d+)')

    RE_FACET_RECORD = re.compile( '\\n *(?P<serial_1>\d+) +(?P<atom_1>\w+) +(?P<res_1>\w+) +(?P<chain_1>[A-Z]?) *(?P<resnr_1>\d+)'+\
                                  ' +(?P<serial_2>\d+) +(?P<atom_2>\w+) +(?P<res_2>\w+) +(?P<chain_2>[A-Z]?) *(?P<resnr_2>\d+)' +\
                                  ' +(?P<pairtype>\w+)' +\
                                  '*\\n(?:[0-9\.\-]+ +){6} (?P<area>[0-9\.]+)')

    
    def __init__( self, model, cr=[0], cl=None, mode=2, breaks=0,
                  catch_err=1, **kw ):
        """
        Create a new Intervor instance for a given protein-protein complex.
        
        @param model: Structure of receptor, ligand and water
        @type  model: Biskit.PDBModel
        @param cr: receptor chains (default: [0] = first chain)
        @type  cr: [ int ]
        @param cl: ligand chains (default: None = all remaining protein chains)
        @type  cl: [ int ]
        @param breaks: consider chain breaks (backbone gaps) (default: 0)
        @type  breaks: bool or 1|0
        @param mode: what to calculate (default 2, = all with shelling order)
        @type  mode: int
        @param catch_err: deviate STDERR to temporary file (default 1)
        @type  catch_err: bool or 0|1
        @param **kw: any other keyword=value pair recognized by Executor
        """
        Executor.__init__( self, 'intervor', catch_err=catch_err, **kw )
        
        assert isinstance( model, B.PDBModel ), \
               'requires PDBModel instance'
        assert model is not None, 'requires PDBModel instance'
        
        self.model = model
        self.breaks= breaks
        self.chains_rec = cr
        self.chains_lig = cl or self.__getLigandChains( model, cr )
        
        self.local_model = None   #: will hold modified copy of model
        
        self.mode = mode
        
        ## intervor puts several output files into current working directory
        ## but respect cwd from Executor.__init__ or ExeConfig/exe_intervor.dat
        self.cwd = self.cwd or tempfile.gettempdir()

        #: will be used by intervor for different output files
        self.f_prefix = tempfile.mktemp( dir=self.exe.cwd )
        self.f_pdb = self.f_prefix + '_intervor.pdb'
        
        self.result = {}
        
 
    def __getLigandChains( self, model, rec_chains ):
        """
        Determine which chains of a model belong to the ligand.
        @return: chain numbers of ligand
        @rtype : [ int ]
        """
        assert isinstance( model, B.PDBModel )

        chains = range( model.lenChains( breaks=self.breaks ) )
        
        ## mark every chain containing at least 1 protein residue with 1
        chain_mask_prot = model.atom2chainMask( model.maskProtein(),
                                            breaks=self.breaks )
        
        r = [ c for c in chains\
              if chain_mask_prot[c] and not c in rec_chains ]
        
        assert len(r) > 0, "Couldn't find any non-receptor protein chains"
        
        return r
    

    def __chainIds( self, model ):
        """
        @return: chain ids of receptor and ligand for intervor, eg. ('AB', 'C')
        @rtype : ( str, str )
        """
        ## get numpy array with chain_id of first atom of each chain
        chainids = N.take( model['chain_id'], 
                           model.chainIndex( breaks=self.breaks ), 0 )
        
        ## extract receptor and ligand chain ids and concat them into string
        ids_rec = ''.join( N.take( chainids, self.chains_rec, 0 ) )
        ids_lig = ''.join( N.take( chainids, self.chains_lig, 0 ) )
        
        return ids_rec, ids_lig


    def __renameWaters( self, model ):
        """Intervor doesn't recognize TIP3 water entries -- rename to HOH"""

        i_h2o = N.flatnonzero( model.maskH2O() )  ## water atoms (indices)
        
        for i in i_h2o:
            model['residue_name'][i] = 'HOH'


##     def __choinWaterChains( self, m ):
##      """workaround the fact that each water gets its own chain"""
##      m._chainIndex = m.chainIndex( breaks=self.breaks )
##      print 'DEBUG', self.chains_lig
##      m._chainIndex = m._chainIndex[ : self.chains_lig[-1]+2 ]


    def prepare( self ):
        """
        Prepare PDB input file for Intervor. Called before program execution.
        * label chains consequtively from 'A' to ..
        * rename TIP3 or any other recognized water into 'HOH'
        * renumber all atoms consequtively starting with 0
        """
        m = self.model.clone()

        self.__renameWaters( m )

        m.addChainId()
        rec_ids, lig_ids = self.__chainIds( m )

        ## last intervor version rejects non-1.0 Occupancies
        m['occupancy'] = N.ones( len(m) ) 
        
        ## intervor reports serial number in facets --> make equal to index
        ## m['serial_number'] = range( len( m ) )  ## use option -p

        ## override commandline args that will be used by run()
        self.args = "-f %s -o %i -C %s -C %s -v 2 -D ." %\
            (self.f_pdb, self.mode, rec_ids, lig_ids) 
        self.args = "-f %s -C %s -C %s -D -p." %\
            (self.f_pdb, rec_ids, lig_ids) 
        
        m.writePdb( self.f_pdb )


    def cleanup( self ):
        """
        Tidy up the mess we created. Called after program execution.
        """
        Executor.cleanup( self )
        
        if not self.debug:
            ## remove all files created by intervor
            for f in glob.glob( self.f_prefix + '*' ):
                T.tryRemove( f )
    

    def isFailed( self ):
        """
        Overrides Executor method
        """
        return self.returncode != 1


    def __parseSO( self, fname ):
        """
        Parse intervor shelling order output file.
        
        @param fname: output file (*_intervor_SO.txt)
        @type  fname: str
        @return: raw table (list of dicts) with intervor output
        @rtype: [ {'edge_index':int, at1_serial:int, ... } ]
        """
        result = []

        try:
            try:
                l = f = None
                f = open( fname )

                ## handle file header (extract value keys from header line)
                l = f.readline()
                while not self.RE_SECTION_SO.match(l):
                        l = f.readline()

                ## handle table (parse values into one dictionary for each line)
                while f:
                    l = f.readline()
                    if l == 'Atom serial#, MinSO, MaxSO\n':
                        break

                    result.append( int( l.split()[1]) )

            except IOError, why:
                raise IntervorError( 'Error accessing intervor output file %r: %r'\
                      % (fname, why))
            except Exception, why:
                raise IntervorError, 'Error parsing intervor output file %r: %r'\
                      % (fname, why) + '\nThe offending line is: ' + repr(l)

        finally:
            try:
                f.close()
            except:
                pass

        return result


    def __parseFacetRecord( self, r ):
        r = self.RE_FACET_RECORD.match( r )
        try:
            d = r.groupdict()

            for k in d:
                if 'serial' in k or 'resnr' in k:
                    d[k] = int( d[k] )

            d['area'] = round( float( d['area'] ), 2 )

            return d
            
        except AttributeError:
            raise IntervorError, 'couldnt parse facet record %r' % r

    def __parseFacets( self, fname ):
        """
        -> Biskit.DictList
        """
        result = []
        
        try:
            r = f = None
            f = open( fname )

            lines = f.readlines()
            lines = ''.join( lines )

            records = lines.split( 'FBIeg' )
            records = [ r for r in records if len( r ) > 10 ]

            result = B.DictList( [ self.__parseFacetRecord(r) for r in records ] )

        except IOError, why:
            raise IntervorError( 'Error accessing intervor output file %r: %r'\
                      % (fname, why))
        except Exception, why:
            raise IntervorError, 'Error parsing intervor output file %r: %r'\
                  % (fname, why) + '\nThe offending line is: ' + repr(r) + '\n'

        return result
    

    def __mapfacets2Atoms( self, raw_table ):
        """
        Create an atom profile 'facets' that contains the facet
        dictionaries each atom is involved in ( empty list if the atom is not
        part of the interface).
        The entries of the profile are Biskit.DictList objects that function
        like normal lists but have additional functions to directly extract 
        values of the dictionaries in the list.
        """
        facet_profile = [ B.DictList() for i in self.model.atomRange() ]
        
        for facet in raw_table.iterDicts():
            a1 = facet['serial_1']
            a2 = facet['serial_2']
            

            facet_profile[a1] += [ facet ]
            facet_profile[a2] += [ facet ]
        
        self.model.atoms.set('facets', facet_profile, comment=\
                             'intervor facets this atom is involved in' )

        self.model['n_facets'] = [ len(f or []) for f in facet_profile ]
        self.model['n_facets','comment'] = 'number of interface facets'

        self.model['so_min'] = [ min( f.valuesOf('so') or [0] )
                                 for f in facet_profile ]
        self.model['so_max'] = [ max( f.valuesOf('so') or [0] )
                                 for f in facet_profile ]
        self.model['facet_area'] = [ N.sum( f.valuesOf('area') or [0])
                                     for f in facet_profile ] 

        self.model['so_min','comment'] = 'smallest Intervor shelling order of '\
                                         + 'any facet this atom is involved in'

        self.model['so_max','comment'] = 'highest Intervor shelling order of '\
                                         +'any facet this atom is involved in'

        self.model['facet_area','comment']='sum of area of all intervor facets'\
                                            + ' this atom is involved in'
        

##     def __createProfiles( self, model ):
##         """
##         Create shelling order profiles from raw facet profile.
##         """
        
##         assert 'facets' in model.atoms, "profile 'facets' not found"
        
##         so_min = [ min( f.valuesOf( 'so' ) or [0] ) for f in model['facets'] ]
##         so_max = [ max( f.valuesOf( 'so' ) or [0] ) for f in model['facets'] ]

##         facet_area = [ N.sum( f.valuesOf( 'facet_area') or [0])
##                        for f in model['facets'] ]
        
##         model['so_min'] = so_min
##         model['so_max'] = so_max
##         model['facet_area'] = facet_area
        


    def finish( self ):
        """Called after a successful intervor run. Overrides Executor hook.
        """
        Executor.finish( self )
        f_so = self.f_prefix + '_intervor_SO.txt'
        f_fac= self.f_prefix + '_intervor_interface-IV.txt'

        facets = B.ProfileCollection()

        facets['so'] = N.array( self.__parseSO( f_so ), N.int32 )

        d = self.__parseFacets( f_fac )

        for k in d.keys():
            facets[ k ] = d.valuesOf( k )
            assert len(facets[k]) == len(facets['so']),\
                   'number of facet records and shelling order records dont match'

        self.result = facets

        self.__mapfacets2Atoms( facets )
##         self.__createProfiles( self.model )
    
        
    def visualize( self, profile='so_min', wat=None, xwat=None ):
        """
        Create a receptor / ligand / water view with a profile color-coded
        on it.
        
        @param profile: name of the profile to color-project ['so_min']
        @type  profile: str
        @param wat: chain ids of water molecules (default: None)
                    These will end up as spheres in group 'water'
        @type  wat: [ int ]
        @param xwat: chain ids of crystal water molecules (default: None)
                     These will end up as spheres in group 'xwater'
        @type  xwat: [ int ]

        @return: Pymoler instance with receptor and ligand colored by shelling
                 order (use pm.run() to make it show up)
        @rtype: Pymoler
        """
        import Biskit.Pymoler as Pymoler
        model = self.model.clone()
        
        ## replace all 0 values by -1 to distinguish them better by color
        model[profile] += (model[profile] == 0) * -1 
        
        rec = model.takeChains( self.chains_rec, breaks=self.breaks )
        lig = model.takeChains( self.chains_lig, breaks=self.breaks )

        water = xwater = None
        
        if wat:
            water  = model.takeChains(wat, breaks=self.breaks)
            water = water.compress( N.greater( water['n_facets'], 0) )

        if xwat:
            xwater = model.takeChains( xwat, breaks=self.breaks ) ## X-ray determined water
            ## kick out waters that do not belong to the interface
            xwater = xwater.compress( N.greater( xwater['n_facets'], 0) )
        
        pm = Pymoler()
        
        pm.addPdb( model, 'com' )
        pm.addPdb( rec, 'rec' )
        pm.addPdb( lig, 'lig' )

        pm.colorAtoms( 'rec', rec[profile])
        pm.colorAtoms( 'lig', lig[profile])

        if xwater:
            pm.addPdb( xwater, 'x-water' )
            pm.colorAtoms( 'x-water', xwater[profile] )

        if water:
            pm.addPdb( water, 'water' )
            pm.colorAtoms( 'water', water[profile] )
        
        pm.add( 'hide' )
        pm.add( 'show surface, rec' )
        pm.add( 'show mesh, lig' )
        pm.add( 'show spheres, x-water' )
        pm.add( 'show spheres, water' )
        pm.add( 'zoom' )

        return pm

        
######## TESTING #########
## Proper test case wrapping comes later
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test class"""
    
    TAGS = [BT.EXE]

    def prepare( self ):
        self.m = B.PDBModel( T.testRoot() + '/com/1BGS.pdb' )

    def test_Intervor( self ):
        """Intervor parsing test"""
        self.x = Intervor( self.m, [0], debug=self.DEBUG, verbose=self.local,
                           catch_err=1 )
        self.x.run()

        if self.local:
           pm = self.x.visualize( profile='so_min', xwat=[2], wat=[3] )
           pm.run()

class TestDry(BT.BiskitTest):
    
    def prepare( self ):
        self.m = B.PDBModel( T.testRoot() + '/com/1BGS.pdb' )

    def test_IntervorDry(self):
        self.x = Intervor( self.m, [0], debug=self.DEBUG, verbose=self.local,
                           catch_err=1, strict=0, validate=0 )
        
        self.x.f_prefix = T.testRoot() + '/intervor/1vfb'

        self.x.finish()

if __name__ == '__main__':

    BT.localTest(debug=0)

##     pm = x.visualize( profile='so_min', xwat=[2], wat=[3] )
##     pm.run()



## new usage:: intervor.exe -f pdbfile -o 2 -C AB -C C -v 2 -D .
