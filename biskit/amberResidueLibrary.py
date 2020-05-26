##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2018 Raik Gruenberg
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

"""
Collect and index AmberResidueType instances from amber topology files.
"""

from biskit import PDBModel, AmberPrepParser, StdLog
import biskit.tools as T

class AmberResidueLibraryError( Exception ):
    pass

class AmberResidueLibrary:
    """
    A collection of reference residue types taken from Amber topology files.
    By default, the collection is initialized from the four all-atom 
    library or "off files" in Biskit/data/amber/residues. 
    The idea here is that the residues of some PDB structure can be matched 
    against this library in order to assign charges or other properties.
    
    Matching residue types are indentified by their 'atom key', which is simply
    the concatenation of a residue's atom names in alphabetical order. This
    means residues are matched by atom content, not by residue name. That's
    important for modified residues, for example, a C-terminal ALA (with an
    additional oxygen) will be matched to the AmberResidueType CALA rather than 
    to ALA. 

    The default all-atom topologies include hydrogen atoms. Structures without
    hydrogens will *not* match. You can add hydrogens with the Reduce class
    (biskit.reduce). 

    Atom names need to conform to Amber conventions -- 
    this can be ensured with `PDBModel.xplor2amber()`.
 
    Use
    ===
    
    >>> model = PDBModel('mystructure.pdb')
    >>> residue = model.resModels()[0]  ## fetch first residue
    >>>
    >>> lib = AmberResidueLibrary()
    >>> refres = lib[ residue ]  ## look up matching AmberResidueType
    >>>
    >>> ## or alternatively:
    >>> refres = lib['all_amino03', 'ALA']
    
    """
    ## list of Amber topology files in decending priority
    F_RESTYPES = ['all_amino03.in',
                  'all_aminoct03.in',
                  'all_aminont03.in',
                  'all_nuc02.in' ]
    

    def __init__(self, topofiles=F_RESTYPES, log=None, verbose=False):
        """
        :param topofiles: list of topology file names \
                          (default `all_*in` in `Biskit/data/amber/residues`)
        :type  topofiles: [ str ]
        :param log: optional LogFile instance (default STDOUT)
        :type  log: biskit.LogFile
        :param verbose: add messages to log (default False)
        :type  verbose: bool
        """       
        self.aindex = {}  ## residue types indexed by atom key
        self.topoindex = {} ## residue types indexed by topo and resname
        
        self.log = log or StdLog()
        self.verbose = verbose
        
        for f in topofiles:
            self.addTopology( f )
        

    def addTopology(self, topofile, override=False):
        """
        Include an additional topology (off) library in the collection.
        
        :param topofile: file name of topology, either full path or
                         simple file name which will then be looked for in 
                         biskit/data/amber/residues.
        :type  topofile: str
        
        :param override: override topologies or residue entries with same name
                         (default False)
        :type  override: False
        
        :return: dictionary of all residue types parsed from topofile indexed
                 by three-letter residue name
        :rtype : {str : AmberResidueType}
        
        :raise: AmberResidueLibraryError if override==False and a topology or
                a residue with identical atom content have already been 
                registered.
        """
        fbase = T.stripFilename( topofile )

        if fbase in self.topoindex and not override:
            raise AmberResidueLibraryError('duplicate topology '+fbase)

        if self.verbose:
            self.log.add('parsing %s...' % topofile )

        resdic = AmberPrepParser( topofile ).residueDict()
        
        if self.verbose:
            self.log.add( 'Read %i residue definitions.\n' % len(resdic) )
        
        self.topoindex[ fbase ] = resdic
        
        for resname, restype in resdic.items():
            akey = restype.atomkey(compress=False)
            
            if akey in self.aindex and not override:
                raise AmberResidueLibraryError('duplicate residue entry: %s -> %s' %\
                      (resname, self.aindex[akey].code))

            self.aindex[ akey ] = restype
        
        return self.topoindex[ fbase ]


    def atomkey( self, residue ):
        """
        Create a string key encoding the atom content of residue.
        
        :param residue: model or AmberResidue
        :type  residue: PDBModel or AmberResidue
        :return: key formed from alphabetically sorted atom content of residue
        :rtype: str
        """
        return residue.atomkey(compress=False)
  
    
    def byAtoms(self, akey, default=None ):
        """
        Identify a matching reference residue by atom content.
        
        :param akey: atomkey or PDBModel with a single residue
        :type  akey: str or PDBModel

        :return: matching reference residue OR None
        :rtype: AmberResidueType
        """
        if isinstance( akey, PDBModel ):
            akey = akey.atomkey(compress=False)
        return self.aindex.get(akey, default)
    
    
    def byName(self, rescode, topo=None ):
        """
        Identify matching reference residue by residue name. Note: residue
        names are not guaranteed to be unique if several topology files have
        been read in (the default set of Amber topologies uses unique names
        though). The optional topo parameter can be used to specify in 
        which topology the residue is looked up.
        
        Note: residue 3 letter names are all UPPERCASE.
        
        :param rescode: three-letter name of residue to look up
        :type  rescode: str
        :param topo: optional (file) name of topology (see also: `topokeys()` )
        :type  topo: str

        :return: matching reference residue
        :rtype: AmberResidueType
        
        :raise: KeyError if the topology or residue name are not found
        """
        if topo:
            fbase = T.stripFilename( topo )
            return self.topoindex[ fbase ][ rescode ]

        for topo, residues in self.topoindex.items():
            if rescode in residues:
                return residues[rescode]

        raise KeyError('No residue type found for name '+str(rescode))
        
    
    def __len__( self ):
        return len(self.aindex)
    
    def __getitem__( self, key ):
        """
        Examples: 
        
            - `reslib[ PDBModel ]` -> `ResidueType` for matching residue
            - `reslib[ str(atomkey) ]` -> `ResidueType` with same atom key
        """
        if type(key) is str:
            return self.aindex[key]
        
        if isinstance(key, PDBModel):
            return self.aindex[ key.atomkey(compress=False) ]
 
        
    def topokeys( self ):
        return list(self.topoindex.keys())
    
    def keys( self ):
        return list(self.aindex.keys())
    
    def values( self ):
        return list(self.aindex.values())
    
#############
##  TESTING        
#############
import biskit.test as BT

class Test(BT.BiskitTest):
    """Test class"""

    def test_amberResidueLibrary( self ):
        """AmberResidueLibrary test"""
        self.lib = AmberResidueLibrary(verbose=self.local)
        ala = self.lib.byName('ALA', 'all_amino03')
        r   = self.lib[ ala ]
        self.assertEqual( r, ala )
        r   = self.lib.byName( 'ALA' )
        self.assertEqual( r, ala )
        
        self.assertEqual( len( self.lib ), 114 )
        
if __name__ == '__main__':

    BT.localTest(debug=False)
