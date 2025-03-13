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
Read Amber residue topology from Amber "off" libraries.
"""
import io
import numpy as N
import os.path as osp
import re

import biskit.tools as T
import biskit as B
import biskit.molUtils as M
 
class AmberPrepError( Exception ):
    pass

class AmberPrepChargeMissmatch( AmberPrepError ):
    pass

class AmberResidueType( B.PDBModel ):
    """
    Standard description of a certain class of residues.
    
    In addition to the normal PDBModel things, this class defines three
    additional fields:
        * name   ... the full long residue name (str), e.g. 'alanine'
        * code   ... the three-letter name used in PDB files (str), e.g. 'ALA'
        * letter ... the single-letter residue code (str), e.g. 'A'
        
    The order and names of atoms are supposed to serve as a reference to check
    and normalize actual residues parsed from PDB files.
    
    Currently, AmberResidueTypes can only be created from Amber prep library
    files.
    
    An alternative method would be to create AmberResidueTypes from some 
    standard PDB file::
    
        >>> m = PDBModel('standard_aa.pdb')
        >>> standard_res = [ AmberResidueType( res ) for res in m.resModels() ]
    
    .. seealso:: `biskit.AmberPrepParser`
    """
    
    def __init__(self, name=None, code=None, letter=None, source=None  ):
        """
        :param name: full residue name (converted to lower letters)
        :type  name: str
        :param code: three-letter residue code (converted to upper letters)
        :type  code: str
        :param letter: single-letter residue code (converted to upper letter)
        :type  letter: str
        :param source: pdb file, model or structure of residue from which to 
                       extract all the other data
        :type  source: PDBModel or str
        """
        
        self.name = name      ## 'alanine'
        self.code = code      ## 'ALA'
        self.letter = letter  ## 'A'
        if self.name:
            self.name = self.name.lower()
        if self.code:
            self.code = self.code.upper()
        if self.letter:
            self.letter = self.letter.upper()
        
        B.PDBModel.__init__( self, source=source )
        
        if source:
            self.__fromPDB()
        

    def __fromPDB( self ):
        self.code = self.atoms['name'][0]        
        self.letter = M.singleAA( [ self.code ] )
        
        
    def __str__( self ):
        return '[%s %3s %3i atoms: %s ]' % \
               (self.__class__.__name__, self.code, len(self), self.name )
    
    def __repr__( self ):
        return str( self )
    
    def take( self, i, rindex=None, cindex=None,
          *initArgs, **initKw ):
        """Overriding PDBModel.take to rescue and copy additional fields"""
        r = B.PDBModel.take(self, i, rindex, cindex, *initArgs, **initKw )
        r.name = self.name
        r.code = self.code
        r.letter = self.letter
        return r

       
class AmberPrepParser( object ):
    """
    Parse Amber Residue libraries (off or prep files) which are usually found 
    in amber/dat/leap/prep.
    
    Usage::
    
        p = AmberOffParser( 'all_amino03.in' )
        residues = [ r for r in p.residueTypes() ]
        
    Returns AmberResidue instances that can be handled like a PDBModel
    although the xyz array contains relative Z-coordinates not real cartesian
    ones. Atom profiles of the AmberResidue instance include:
        * charge ... atomic partial charge (float)
        * name   ... atom name (str)
        * amber_type ... Amber atom type (str)
    """
    
    F_DEFAULT = 'all_amino03.in'
    
    def __init__(self, f_in=None ):
        """
        :param f_in: amber "off" or "prep" file with residue definitions
                     if not existing, we will look for a file with this name
                     in the Biskit data folder (data/amber/residues)
                     (default: 'all_amino03.in')
        :type  f_in: str
        """
        f_in = f_in or self.F_DEFAULT
        if not osp.exists( T.absfile( f_in ) ):
            f_in = T.dataRoot() + '/amber/residues/' + f_in
        
        self.firstrecord = re.compile(r'db[0-9]+\.dat' )
        
        self.s = open( T.absfile( f_in ), 'r' ).read()
        self.s = self.firstrecord.split( self.s )[-1] #skip until first residue
  

    def residueBlocks( self ):
        i_from = 0
        i_to = self.s.find( 'DONE', i_from )

        while i_to != -1:
            yield self.s[i_from : i_to]
            i_from = i_to + 4
            i_to   = self.s.find( 'DONE', i_from )
            
    
    def atomLines( self, s ):
        """:param s: str, atom block
        """
        start = s.find('0.0')
        end = s.find( '\n\n', start + 10 )

        h = io.StringIO( s[start:end] )
        for i in range(4):
            h.readline()
        
        line = h.readline().strip()
        while line:
            yield line
            line = h.readline().strip()


    def parseCharges( self, s ):
        s = s[ s.find('CHARGE') : ]
        s = s[ s.find('\n') : ]
        
        endpos1 = s.find('LOOP')
        endpos2 = s.find('IMPROPER')
        endpos = endpos2
        if endpos1 != -1 and endpos1 < endpos2:
            endpos = endpos1
        
        s = s[ : endpos ]
        s = s.strip()
        return s.split()

    
    def parseResidue( self, s ):
        """:param s: str, residue record
        """
        r = {}
        h = io.StringIO( s )
        line = h.readline().strip()
        while not line:
            line = h.readline().strip()

        r['name'] = line
        h.readline()
        r['code'] = h.readline().split()[0]
        
        atoms = {'serial_number': [], 'name':[], 'amber_type':[],
                 'xyz':[], 'partial_charge':[] }
        for l in self.atomLines( s ):
            items = l.split()
            atoms['serial_number'].append( int( items[0] ) )
            atoms['name'].append( items[1] )
            atoms['amber_type'].append( items[2] )
            atoms['xyz'].append( items[7:10] )
            try:
                ## BIG FAT WARNING: charge column sometimes deviates from 
                ## charge record further down -- only use as fall-back
                ## Example mess up: HIP in all_amino03
                atoms['partial_charge'].append( items[10] )
            except:
                pass  ## charge column is not always present
            
        try:
            ## the charge record seems more reliable but is not always there
            q = self.parseCharges( s )
            if len(q) != len( atoms['xyz'] ):
                if r['name'] == 'HISTIDINE PLUS':
                    pass
                raise AmberPrepChargeMissmatch('broken charge record in '+\
                      r['code'])
            atoms['partial_charge'] = q
        except IndexError:
            pass
        except AmberPrepChargeMissmatch as why:
            pass

        if atoms['partial_charge'] == []:
            raise AmberPrepError('failed to parse charges for '+r['code'])

        atoms['partial_charge'] = N.array( atoms['partial_charge'], float )
        atoms['xyz']    = N.array( atoms['xyz'], float )
        return r, atoms

    
    def createResidue( self, res, atoms ):
        """
        """
        r = AmberResidueType( **res )
        r.letter = M.singleAA( [r.code] )[0]
        r.xyz = atoms['xyz']

        n = len(r.xyz)
        for key in B.PDBModel.PDB_KEYS:
            r[key] = n * ['']
        
        r['type'] = ['ATOM'] * n
        r['residue_name'] = [r.code] * n
        r['residue_number'] = [1] * n
        r['occupancy'] = [1.0] * n
        r['after_ter'] = [0] * n
        r['temperature_factor'] = [0] * n
        
        del atoms['xyz']
        for key, profile in atoms.items():
            r.atoms[key] = profile

        return r

    
    def residueTypes( self ):
        for resblock in self.residueBlocks():
            r, atoms = self.parseResidue( resblock )
            yield self.createResidue( r, atoms )
        
    def residueDict( self ):
        """
        :return: dict, residue types indexed by 3-letter residue code,
                 example: {'ala':[AmberResidueType ALA], 'cys':[...}
        """
        r = dict( [(res.code,res) for res in self.residueTypes()] )
        return r
            

#############
##  TESTING        
#############
import biskit.test as BT
import glob

class Test(BT.BiskitTest):
    """Test class"""

    def test_amberPrepParser( self ):
        """AmberPrepParser test"""

        files = glob.glob( T.dataRoot()+'/amber/residues/*in')
        files = [ osp.basename( f ) for f in files ]
        results = {}
        if self.local:
            print()
        
        for f in files:
            if self.local:
                print('working on ', f)
            self.p = AmberPrepParser( f )
            self.r = [ r for r in self.p.residueTypes() ]
            self.assertTrue( len(self.r) > 10 )
        
            if self.local:
                print('\tparsed %i residue types from %s' % (len(self.r), f))
            results[ f ] = self.r
        
        self.assertEqual( len(results['all_amino03.in']), 33 )
        
        dic = self.p.residueDict()
        self.assertTrue( isinstance( list(dic.values())[0], AmberResidueType) )
        
        self.results = results
        if self.local:
            for res in results['all_nuc02.in']:
                print(res)
                

if __name__ == '__main__':

    BT.localTest(debug=False)

    
