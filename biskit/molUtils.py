## numpy-oldnumeric calls replaced by custom script; 09/06/2016
## Automatically adapted for numpy-oldnumeric Mar 26, 2007 by alter_code1.py

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
Utilities for handling structures and sequences
"""
## see: https://www.python.org/dev/peps/pep-0366/
## allow relative imports when calling module as main script for testing
if __name__ == "__main__" and __package__ is None:
    import biskit
    __package__ = "biskit"

from biskit import EHandler
from biskit import tools as t
from biskit.core import oldnumeric as N0

import copy
import types

class MolUtilError( Exception ):
    pass

#: translate PDB amino acid names to single letter code
aaDicStandard =\
              {'asp':'D', 'glu':'E', 'lys':'K', 'his':'H', 'arg':'R',
               'gln':'Q', 'asn':'N', 'ser':'S', 'asx':'B', 'glx':'Z',
               'phe':'F', 'trp':'W', 'tyr':'Y',
               'gly':'G', 'ala':'A', 'ile':'I', 'leu':'L', 'cys':'C',
               'met':'M', 'thr':'T', 'val':'V', 'pro':'P' }

#: same for nucleic acids (incomplete)
nsDicStandard = {'a':'a', 'g':'g', 'c':'c', 't':'t', 'u':'u',
                 'a3':'a', 'g3':'g', 'c3':'c', 't3':'t', 'u3':'u',
                 'a5':'a', 'g5':'g', 'c5':'c', 't5':'t', 'u5':'u',
                 'da':'a','dg':'g','dc':'c','dt':'t',
                 'da3':'a','dg3':'g','dc3':'c','dt3':'t',
                 'da5':'a','dg5':'g','dc5':'c','dt5':'t'
                 }

#: extend aaDicStandard with non-standard residues
aaDic = copy.copy( aaDicStandard )

aaDic.update( {'cyx':'C', 'hid':'H', 'hie':'H', 'hip':'H',
               'unk':'X', 'ace':'X', 'nme':'X'} )#, 'ndp':'X' } )

#: extend nsDicStandard with non-standard residues
nsDic = copy.copy( nsDicStandard )

nsDic.update( {'atp':'a', 'gtp':'g', 'ctp':'c', 'ttp':'t', 'utp':'u',
               'adp':'a', 'gdp':'g', 'cdp':'c', 'tdp':'t', 'udp':'u',
               'amp':'a', 'gmp':'g',
               'fad':'f', 'fmp':'f',
               'nad':'n',
           } )

#: translate common hetero residues to pseudo single letter code
xxDic = {'tip3':'~', 'hoh':'~', 'wat':'~', 'cl-':'-', 'na+':'+', 'ca':'+',
         'ndp':'X', 'nap':'X'}

#: translate standard PDB amino and nucleic acid names to single letter code
resDicStandard = copy.copy( aaDicStandard )
resDicStandard.update( nsDicStandard )

#: extend resDicStandard with common non-standard names
resDic = copy.copy( aaDic )
resDic.update( nsDic )
resDic.update( xxDic )

## map non-standard amino acid names to closest standard amino acid
##
## Data from: http://www.ccp4.ac.uk/html/lib_list.html#peptide_synonyms
## More info at: http://xray.bmc.uu.se/hicup/XXX/ where XXX is the residue code
##
## NOT ADDED:
##     SAR  SARCOSINE
##     PCA  5-pyrrolidone-2-carboxylic_acid
##     INI  Amidinated_lysine_with_methyl_isonicotinimida
##     SAH  S-ADENOSYL-L-HOMOCYSTEINE
##     SAM  S-ADENOSYLMETHIONINE
##     LLP  LYSINE-PYRIDOXAL-5*-PHOSPHATE

##     ACE  acetyl 
##     FOR  Formyl
##     BOC  TERT-BUTYLOXYCARBONYL GROUP
##     MLE  N-METHYLLEUCINE
##     MVA  N-METHYLVALINE
##     IVA  Isovaleric_acid  
##     STA  STATINE
##     ETA  ethanolamine
##     TFA  TRIFLUOROACETYL GROUP
##     ANI  4-TRIFLUOROMEHYLANILINE
##     MPR  BETA-MERCAPTOPROPIONATE
##     DAM  N-METHYL-ALPHA-BETA-DEHYDROALANINE
##     ACB  2-AMINO-3-CARBONYLBUTANOIC ACID
##     ADD  2,6,8-TRIMETHYL-3-AMINO-9-BENZYL-9-M
##     CXM  N-CARBOXYMETHIONINE
##     DIP  DIPENTYLAMINE
##     BAL   BETA-ALANINE

nonStandardAA={ 'UNK':'ALA', 'ABA':'ALA', 'B2A':'ALA',
                'ORN':'ARG',
                'ASX':'ASP',
                'CSH':'CYS', 'OCS':'CYS', 'CSO':'CYS',
                'GLX':'GLU', 'CGU':'GLU', 'ILG':'GLU',
                'B2I':'ILE',
                'BLE':'LEU',
                'KCX':'LYS', 'BLY':'LYS',
                'MSE':'MET',
                'B1F':'PHE', 'B2F':'PHE',
                'HYP':'PRO', '5HP':'PRO',
                'SEP':'SER',
                'TYS':'TYR',
                'B2V':'B2V',
                'HIE':'HIS', 'HID':'HIS', 'HIP':'HIS',
                'CYX':'CYS' }

#: heavy atoms of amino acids in standard order, OXT applies to C term only
aaAtoms={'GLY':['N','CA','C','O', 'OXT' ],
         'ALA':['N','CA','C','O', 'CB', 'OXT'],
         'VAL':['N','CA','C','O','CB','CG1','CG2', 'OXT'],
         'LEU':['N','CA','C','O','CB','CG','CD1','CD2', 'OXT'],
         'ILE':['N','CA','C','O','CB','CG1','CG2','CD1', 'OXT'],
         'MET':['N','CA','C','O','CB','CG','SD','CE', 'OXT'],
         'PRO':['N','CA','C','O','CB','CG','CD', 'OXT'],
         'PHE':['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ',
                'OXT'],
         'TRP':['N','CA','C','O','CB','CG','CD1','CD2','NE1','CE2','CE3',
                'CZ2','CZ3','CH2', 'OXT'],
         'SER':['N','CA','C','O','CB','OG', 'OXT'],
         'THR':['N','CA','C','O','CB','OG1','CG2', 'OXT'],
         'ASN':['N','CA','C','O','CB','CG','OD1','ND2', 'OXT'],
         'GLN':['N','CA','C','O','CB','CG','CD','OE1','NE2', 'OXT'],
         'TYR':['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ','OH',
                'OXT'],
         'CYS':['N','CA','C','O','CB','SG', 'OXT'],
         'LYS':['N','CA','C','O','CB','CG','CD','CE','NZ', 'OXT'],
         'ARG':['N','CA','C','O','CB','CG','CD','NE','CZ','NH1','NH2', 'OXT'],
         'HIS':['N','CA','C','O','CB','CG','ND1','CD2','CE1','NE2', 'OXT'],
         'ASP':['N','CA','C','O','CB','CG','OD1','OD2', 'OXT'],
         'GLU':['N','CA','C','O','CB','CG','CD','OE1','OE2', 'OXT']}

#: dictionary of elements
elements = { 'carbon':['C', 'CD2', 'CZ2', 'CB', 'CA', 'CG', 'CE', 'CD', 'CZ',
                       'CH2', 'CE3', 'CD1', 'CE1', 'CZ3', 'CG1', 'CG2', 'CE2'],
             'nitrogen':['NZ', 'ND2', 'NH1', 'NH2', 'ND1', 'NE1', 'NE2',
                         'NE', 'N'],
             'oxygen':['OG', 'OE2', 'OXT', 'OD1', 'OE1', 'OH', 'OG1', 'OD2',
                       'O'],
             'suplphur':['SG', 'SD'],
             'clustering_BDZ':['C','CB','CD','CD1','CD2','CZ','CZ2','CZ3',
                               'ND1','ND2','NZ','OD1','OD2','SD' ],
             'clustering_ABDZ':['C','CA','CB','CD','CD1','CD2','CZ','CZ2',
                                'CZ3',
                               'ND1','ND2','NZ','OD1','OD2','SD' ],
             'clustering_G':['C','CG','CG1','OG','OG1','SG' ],
             'clustering_B':['C','CB'],
             'clustering_AG':['C','CA','CG','CG1','OG','OG1','SG' ],
             'clustering_AGE':['C','CA','CG','CG1','OG','OG1','SG','NE','OE1',
                               'CE1','CE','CE3' ],
             'clustering_BD':['C','CB','CD','CD1','OD1','SD' ],
             'clustering_ABD':['C','CA','CB','CD','CD1','OD1','SD' ],
             'clustering_AB':['C','CA','CB']}

#: number of attached H for each heavy atom in each amino acid
aaAtomsH={'XXX':{'N':1,'CA':1,'C':0,'O':0,'OXT':0},
          'GLY':{},
          'ALA':{'CB':3},
          'VAL':{'CB':0,'CG1':3,'CG2':3},
          'LEU':{'CB':2,'CG':0,'CD1':3,'CD2':3},
          'ILE':{'CB':0,'CG1':1,'CG2':3,'CD1':3},
          'MET':{'CB':2,'CG':2,'SD':0,'CE':3 },
          'PRO':{'N':0,'CB':2,'CG':2,'CD':2},
          'PHE':{'CB':2,'CG':0,'CD1':1,'CD2':1,'CE1':1,'CE2':1,'CZ':1},
          'TRP':{'CB':2,'CG':0,'CD1':1,'CD2':0,'NE1':1,'CE2':0,'CE3':1,
                 'CZ2':1,'CZ3':1,'CH2':1},
          'SER':{'CB':2,'OG':1},
          'THR':{'CB':0,'OG1':1,'CG2':3},
          'ASN':{'CB':2,'CG':0,'OD1':0,'ND2':2},
          'GLN':{'CB':2,'CG':2,'CD':0,'OE1':0,'NE2':2},
          'TYR':{'CB':2,'CG':0,'CD1':1,'CD2':1,'CE1':1,'CE2':1,'CZ':0,'OH':1},
          'CYS':{'CB':2,'SG':1},
          'LYS':{'CB':2,'CG':2,'CD':2,'CE':2,'NZ':3},
          'ARG':{'CB':2,'CG':2,'CD':2,'NE':1,'CZ':0,'NH1':2,'NH2':2},
          'HIS':{'CB':2,'CG':0,'ND1':1,'CD2':1,'CE1':1,'NE2':0},
          'ASP':{'CB':2,'CG':0,'OD1':0,'OD2':0},
          'GLU':{'CB':2,'CG':2,'CD':0,'OE1':0,'OE2':0} }

for aa in aaAtomsH:
    default = copy.copy( aaAtomsH['XXX'] )
    default.update( aaAtomsH[aa] )
    aaAtomsH[aa] = default

## work in progress...heavy atoms of nucleic acids in standard order
nsAtoms={
    'ATP':['PG', 'O1G', 'O2G', 'O3G', 'PB', 'O1B', 'O2B', 'O3B', 'PA', 'O1A',
           'O2A', 'O3A', 'O5*', 'C5*', 'C4*', 'O4*', 'C3*', 'O3*', 'C2*',
           'O2*', 'C1*', 'N9', 'C8', 'N7', 'C5', 'C6', 'N6', 'N1', 'C2', 'N3',
           'C4'],
    'GTP':['PG', 'O1G', 'O2G', 'O3G', 'PB', 'O1B', 'O2B', 'O3B', 'PA', 'O1A',
           'O2A', 'O3A', 'O5*', 'C5*', 'C4*', 'O4*', 'C3*', 'O3*', 'C2*',
           'O2*', 'C1*', 'N9', 'C8', 'N7', 'C5', 'C6', 'O6', 'N1', 'C2', 'N2',
           'N3', 'C4'],
    'DA': ['P', 'O1P', 'O2P', "O5'", "C5'", "H5'1", "H5'2", "C4'", "H4'", "O4'",
           "C1'", "H1'", 'N9', 'C8', 'H8', 'N7', 'C5', 'C6', 'N6', 'H61', 'H62',
           'N1', 'C2', 'H2', 'N3', 'C4', "C3'", "H3'", "C2'", "H2'1", "H2'2", 
           "O3'"],
    'DC': ['P', 'O1P', 'O2P', "O5'", "C5'", "H5'1", "H5'2", "C4'", "H4'", "O4'",
           "C1'", "H1'", 'N1', 'C6', 'H6', 'C5', 'H5', 'C4', 'N4', 'H41', 'H42',
           'N3', 'C2', 'O2', "C3'", "H3'", "C2'", "H2'1", "H2'2", "O3'"],
    'DG': ['P', 'O1P', 'O2P', "O5'", "C5'", "H5'1", "H5'2", "C4'", "H4'", "O4'",
           "C1'", "H1'", 'N9', 'C8', 'H8', 'N7', 'C5', 'C6', 'O6', 'N1', 'H1',
           'C2', 'N2', 'H21', 'H22', 'N3', 'C4', "C3'", "H3'", "C2'", "H2'1",
           "H2'2", "O3'"],
    'DT': ['P', 'O1P', 'O2P', "O5'", "C5'", "H5'1", "H5'2", "C4'", "H4'", "O4'",
           "C1'", "H1'", 'N1', 'C6', 'H6', 'C5', 'C7', 'H71', 'H72', 'H73', 
           'C4', 'O4', 'N3', 'H3', 'C2', 'O2', "C3'", "H3'", "C2'", "H2'1", 
           "H2'2", "O3'"],
    'RA': ['P', 'O1P', 'O2P', "O5'", "C5'", "H5'1", "H5'2", "C4'", "H4'", "O4'",
           "C1'", "H1'", 'N9', 'C8', 'H8', 'N7', 'C5', 'C6', 'N6', 'H61', 'H62',
           'N1', 'C2', 'H2', 'N3', 'C4', "C3'", "H3'", "C2'", "H2'1", "O2'",
           "HO'2", "O3'"],
    'RC': ['P', 'O1P', 'O2P', "O5'", "C5'", "H5'1", "H5'2", "C4'", "H4'", "O4'",
           "C1'", "H1'", 'N1', 'C6', 'H6', 'C5', 'H5', 'C4', 'N4', 'H41', 'H42',
           'N3', 'C2', 'O2', "C3'", "H3'", "C2'", "H2'1", "O2'", "HO'2", "O3'"],
    'RG': ['P', 'O1P', 'O2P', "O5'", "C5'", "H5'1", "H5'2", "C4'", "H4'", "O4'",
           "C1'", "H1'", 'N9', 'C8', 'H8', 'N7', 'C5', 'C6', 'O6', 'N1', 'H1',
           'C2', 'N2', 'H21', 'H22', 'N3', 'C4', "C3'", "H3'", "C2'", "H2'1",
           "O2'", "HO'2", "O3'"],
    'RU': ['P', 'O1P', 'O2P', "O5'", "C5'", "H5'1", "H5'2", "C4'", "H4'", "O4'",
           "C1'", "H1'", 'N1', 'C6', 'H6', 'C5', 'H5', 'C4', 'O4', 'N3', 'H3', 
           'C2', 'O2', "C3'", "H3'", "C2'", "H2'1", "O2'", "HO'2", "O3'"],
    'MG' :['MG'],
    'NDP':['P1', 'O1', 'O2', 'O5R', 'C5R', 'O1R', 'C4R', 'C3R', 'O3R', 'C2R',
           'O2R', 'C1R', 'N9', 'C8', 'N7', 'C5', 'C6', 'N6', 'N1', 'C2',
           'N3', 'C4', 'O10', 'P2', 'O11', 'O21', 'O51R', 'C51R', 'O11R',
           'C41R', 'C31R', 'O31R', 'C21R', 'O21R', 'C11R', 'N11', 'C61',
           'C51', 'C71', 'O71', 'N71', 'C41', 'C31', 'C21', 'P3', 'O3',
           'O4', 'O5', 'H8', 'H9', 'H7', 'H6', 'H1', 'H5', 'H4', 'H13',
           'H11', 'H12', 'H10', 'H18', 'H19', 'H17', 'H16', 'H3', 'H15',
           'H2', 'H14', 'H23', 'H24', 'H25', 'H22', 'H26', 'H21', 'H20'] }

for res in ['DA','DC','DG','DT','RA','RC','RG','RU']:
    #delete H
    nsAtoms[ res ] = [ a for a in nsAtoms[res] if a[0] != 'H' ]

    # create 3' and 5' versions
    nsAtoms[ res + '3' ] = nsAtoms[res] + ['H3T']
    nsAtoms[ res + '5' ] = ['H5T'] + nsAtoms[res][3:]

nsAtoms['NAP'] = nsAtoms['NDP'].remove('H26')

#: map AA and NS and some other residue names to single letter code
resDic = copy.copy( aaDic )
resDic.update( nsDicStandard )

#: map AA and NS residue names to list of allowed heavy atoms
atomDic = copy.copy( aaAtoms )
atomDic.update( nsAtoms )

#: some common synonyms of atom names
atomSynonyms = { "O'":'O', 'OT1':'O', "O''":'OXT', 'OT2':'OXT',
                 'O1':'O', 'O2':'OXT',
                 'CD':'CD1'}

hydrogenSynonyms = { 'H':'HN',      '1HE2':'HE21', '2HE2':'HE22',
                     '1HH1':'HH11', '2HH1':'HH12', '1HH2':'HH21',
                     '2HH2':'HH22', '1HD2':'HD21', '2HD2':'HD22' }

###################
## Hydrogen bond
##
hbonds={ 'donors': {'GLY':['H','H1','H2','H3'],
                    'ALA':['H','H1','H2','H3'],
                    'VAL':['H','H1','H2','H3'],
                    'LEU':['H','H1','H2','H3'],
                    'ILE':['H','H1','H2','H3'],
                    'MET':['H','H1','H2','H3'],
                    'PRO':['H','H1','H2','H3'],
                    'PHE':['H','H1','H2','H3'],
                    'TRP':['H','H1','H2','H3','HE1'],
                    'SER':['H','H1','H2','H3','HG'],
                    'THR':['H','H1','H2','H3','HG1'],
                    'ASN':['H','H1','H2','H3','1HD2','2HD2'],
                    'GLN':['H','H1','H2','H3','1HE2','2HE2'],
                    'TYR':['H','H1','H2','H3','HH'],
                    'CYS':['H','H1','H2','H3','HG'],
                    'LYS':['H','H1','H2','H3','HZ1','HZ2','HZ3'],
                    'ARG':['H','H1','H2','H3','HE','1HH1','2HH1',
                           '1HH2','2HH2'],
                    'HIS':['H','H1','H2','H3','HD1','HE2'],
                    'ASP':['H','H1','H2','H3'],
                    'GLU':['H','H1','H2','H3']},
         'acceptors': {'GLY':['O','OXT' ],
                       'ALA':['O','OXT'],
                       'VAL':['O','OXT'],
                       'LEU':['O','OXT'],
                       'ILE':['O','OXT'],
                       'MET':['O','SD','OXT'],
                       'PRO':['O','OXT'],
                       'PHE':['O','OXT'],
                       'TRP':['O','OXT'],
                       'SER':['O','OG', 'OXT'],
                       'THR':['O','OG1','CG2', 'OXT'],
                       'ASN':['O','OD1','OXT'],
                       'GLN':['O','OE1','OXT'],
                       'TYR':['O','OH','OXT'],
                       'CYS':['O','SG','OXT'],
                       'LYS':['O','OXT'],
                       'ARG':['O','OXT'],
                       'HIS':['O','OXT'],
                       'ASP':['O','OD1','OD2', 'OXT'],
                       'GLU':['O','OE1','OE2', 'OXT']} }


##############################
## Polar hydrogen connectivity -- PARAM19

polarH = {'GLY':{'H':'N','H1':'N','H2':'N','H3':'N'},
          'ALA':{'H':'N','H1':'N','H2':'N','H3':'N'},
          'VAL':{'H':'N','H1':'N','H2':'N','H3':'N'},
          'LEU':{'H':'N','H1':'N','H2':'N','H3':'N'},
          'ILE':{'H':'N','H1':'N','H2':'N','H3':'N'},   
          'MET':{'H':'N','H1':'N','H2':'N','H3':'N'},
          'PRO':{'H':'N','H1':'N','H2':'N','H3':'N'},
          'PHE':{'H':'N','H1':'N','H2':'N','H3':'N'},
          'TRP':{'H':'N','H1':'N','H2':'N','H3':'N',
                 'HE1':'NE1'},
          'SER':{'H':'N','H1':'N','H2':'N','H3':'N',
                 'HG':'OG'},
          'THR':{'H':'N','H1':'N','H2':'N','H3':'N',
                 'HG1':'OG1'},
          'ASN':{'H':'N','H1':'N','H2':'N','H3':'N',
                 'HD21':'ND2','HD22':'ND2'},
          'GLN':{'H':'N','H1':'N','H2':'N','H3':'N',
                 'HE21':'NE2','HE22':'NE2'},
          'TYR':{'H':'N','H1':'N','H2':'N','H3':'N',
                 'HH':'OH'},
          'CYS':{'H':'N','H1':'N','H2':'N','H3':'N'},
          'LYS':{'H':'N','H1':'N','H2':'N','H3':'N',
                 'HZ1':'NZ','HZ2':'NZ','HZ3':'NZ'},
          'ARG':{'H':'N','H1':'N','H2':'N','H3':'N',
                 'HE':'NE', 'HH11':'NH1','HH12':'NH1',
                 'HH21':'NH2','HH22':'NH2'},
          'HIS':{'H':'N','H1':'N','H2':'N','H3':'N',
                 'HD1':'ND1','HE2':'NE2'},
          'ASP':{'H':'N','H1':'N','H2':'N','H3':'N'},
          'GLU':{'H':'N','H1':'N','H2':'N','H3':'N'}}


## Scoring matrix for protein-protein interaction surfaces
## (Volume normalized values, Table IV in reference)
##
## The Matrix is based on data from a db of 621 noneredundant protein-protein
## complexes, a CB-CB (CA for Gly) of 6 A was used
##
## Reference:
## "Residue Frequencies and Pair Preferences at Protein-Protein Interfaces"
##            F. Glaser, D. M. Steinberg, I. A. Vakser and N0. Ben-Tal,
##            Proteins 43:89-102 (2001)
##
## Warning. This is just half of the matrix (above diagonal), the residue names
## in the pairs is sorted in the same order as in Complex.resPairCounts()

pairScore = {'WW': 5.85, 'WY': 6.19, 'RT': 3.77, 'RV': 4.18, 'RW': 8.57, 'RR': 2.87,
             'RS': 2.82, 'RY': 5.28, 'GW': 1.42, 'GV':-0.41, 'GT': 0.21, 'GS':-1.53,
             'GR': 1.59, 'GQ': 1.70, 'GP':-0.51, 'GY': 1.25, 'GG':-4.40, 'GN':-0.54,
             'GM': 0.91, 'GL':-0.37, 'GK': 1.33, 'GI': 0.77, 'GH': 1.08, 'SS':-0.09,
             'IY': 5.61, 'HY': 6.05, 'HR': 4.90, 'HS': 0.80, 'HP': 2.89, 'HQ': 4.00,
             'HV': 3.21, 'HW': 6.46, 'HT': 2.71, 'KN': 3.17, 'HK': 2.72, 'HH': 5.37,
             'HI': 3.38, 'HN': 2.38, 'HL': 4.88, 'HM': 4.65, 'ST': 1.91, 'PR': 3.99,
             'PS': 1.33, 'PP': 0.60, 'PQ': 3.50, 'PV': 2.90, 'PW': 7.87, 'PT': 2.65,
             'PY': 4.22, 'IQ': 3.60, 'IP': 3.27, 'AK': 2.13, 'EM': 3.88, 'EL': 3.12,
             'EN': 2.68, 'EI': 3.20, 'EH': 2.30, 'EK': 5.32, 'EE': 1.65, 'EG':-0.89,
             'EF': 2.87, 'IT': 3.05, 'EY': 4.54, 'ET': 2.88, 'EW': 1.20, 'IV': 4.91,
             'EQ': 1.95, 'EP': 3.17, 'ES': 2.60, 'ER': 5.75, 'II': 3.89, 'MM': 6.02,
             'MN': 2.30, 'AS': 0.39, 'MT': 2.09, 'MW': 4.89, 'MV': 4.37, 'MQ': 4.18,
             'MP': 3.38, 'MS': 1.61, 'MR': 3.62, 'MY': 4.81, 'IL': 4.59, 'FP': 4.25,
             'FQ': 4.25, 'FR': 4.49, 'FS': 1.75, 'FT': 3.34, 'VV': 3.74, 'FV': 4.69,
             'FW': 5.83, 'FY': 5.83, 'AV': 2.57, 'FF': 5.34, 'FG': 0.14, 'FH': 3.47,
             'FI': 5.33, 'FK': 3.57, 'FL': 4.86, 'FM': 5.28, 'FN': 3.11, 'EV': 3.22,
             'NN': 2.92, 'NY': 3.66, 'NP': 3.09, 'NQ': 3.45, 'NR': 3.85, 'NS': 1.77,
             'NT': 2.52, 'NV': 1.36, 'NW': 3.54, 'CK': 2.05, 'CI': 1.76, 'CH': 4.12,
             'CN':-0.42, 'CM': 1.84, 'CL': 2.93, 'CC': 7.65, 'CG':-0.25, 'CF': 3.68,
             'CE': 2.51, 'CD': 0.24, 'CY': 2.47, 'CS': 2.48, 'CR': 2.81, 'CQ': 1.33,
             'CP': 2.47, 'CW': 2.14, 'CV': 2.89, 'CT': 1.03, 'SY': 2.30, 'VW': 2.92,
             'KK': 3.24, 'SW': 2.87, 'SV': 1.42, 'KM': 3.93, 'KL': 3.15, 'KS': 2.74,
             'KR': 2.29, 'KQ': 3.50, 'KP': 3.75, 'KW': 5.76, 'KV': 4.45, 'KT': 3.67,
             'KY': 5.26, 'DN': 3.85, 'DL': 1.40, 'DM': 0.36, 'DK': 3.90, 'DH': 5.20,
             'DI': 2.30, 'DF': 0.99, 'DG':-0.08, 'DD': 0.13, 'DE': 0.08, 'YY': 5.93,
             'DY': 1.76, 'DV': 1.93, 'DW': 2.62, 'DT': 3.88, 'DR': 4.94, 'DS': 2.94,
             'DP': 1.46, 'DQ': 3.26, 'TY': 3.14, 'LN': 2.31, 'TW': 5.12, 'LL': 4.03,
             'LM': 5.32, 'LV': 4.20, 'LW': 5.77, 'LT': 2.07, 'LR': 4.99, 'LS': 1.41,
             'LP': 2.50, 'LQ': 3.46, 'LY': 4.19, 'AA':-0.52, 'AC': 1.46, 'AE': 1.71,
             'AD': 1.13, 'AG':-1.77, 'AF': 3.00, 'AI': 2.84, 'AH': 2.59, 'IS': 1.00,
             'IR': 3.80, 'AM': 2.30, 'AL': 2.77, 'IW': 6.24, 'AN': 1.69, 'AQ': 1.72,
             'AP': 1.22, 'IK': 3.23, 'AR': 1.90, 'IM': 5.25, 'AT': 1.21, 'AW': 3.37,
             'IN': 1.59, 'AY': 2.47, 'VY': 3.95, 'QQ': 2.83, 'QS': 2.00, 'QR': 4.50,
             'QT': 1.82, 'QW': 1.37, 'QV': 3.22, 'QY': 2.05, 'TV': 2.83, 'TT': 1.27}

## various constants
boltzmann    = 1.38066e-23    ## [J/K]
NA           = 6.02214199e+23 ## Avogadro constant [1/mol]
planck2      = 1.0545727e-34  ## [J s], h/2Pi
euler        = N0.e
mu           = 1.66056e-27    ## atomic mass unit in [kg]
angstroem    = 1e-10          ## [m]
calorie      = 4.184          ## [J]

#: dictionary with relative atomic mass of elements {'H':1.01, 'ZN':65.39, ...}
atomMasses = { 'H':1.00797, 'C':12.01115, 'N':14.0067,
               'S':32.064,  'O':15.9994,  'P':30.9738, 'ZN': 65.39 }

def allAACodes():
    """
    :return: list of all single AA codes, including B, Z, X
    :rtype: [str]
    """
    result = []
    for aa in aaDic.values():
        if not aa in result:
            result += aa

    return result


def allAA():
    """
    :return: list of all 20 'exact' single AA codes.
    :rtype: [str]
    """
    result = allAACodes()

    for a in ['Z','B','X']:
        result.remove( a )

    return result


def elementType( eLetter ):
    """
    Classify an atom as polar or unpolar::
      atomType( eLetter ) -> list of types this element belongs to

    :param eLetter: atom name
    :type  eLetter: str
    
    :return: return 'p' for polar, 'u' for unpolar and None if not
             in classified
    :rtype: p|u OR None
    """
    types = {'p' : ['N','O','H','Cl'],  ## polar
             'u' : ['C','S'] }          ## unpolar

    for key, values in types.items():
        if eLetter in values:
            return key
    return None


def resType( resCode ):
    """
    Classify residues as aromatic (a), charged (c) or polar (p).

    :param resCode: amino acid code
    :type  resCode: str
    
    :return: list of types this residue belongs to...
    :rtype: a|c|p OR None
    """
    types = {'a' : ['F','Y','W','H'],     ## aromatic
             'c' : ['E','D','L','R','H'], ## charged
             'p' : ['Q','N','S'] }        ## polar

    result = []

    for t in types.keys():
        if resCode in types[t]:
            result += [t]

    if result == []:
        result = ['u']

    return result


def singleAA(seq, xtable=None, nonstandard=True, unknown='?' ):
    """
    convert list with 3-letter AA code to list with 1-letter code
    
    :param seq: amino acid sequence in 3-letter code
    :type  seq: [str]
    :param xtable: dictionary with additional str:single_char mapping
    :type  xtable: dict
    :param nonstandard: support non-standard residue names (default True)
    :type  nonstandard: bool
    :param unknown: letter to use for unknown residues [default: '?']
    :type unknown: str
    
    :return: list with 1-letter code; C{ ['A','C','L','A'...]}
    :rtype: [str]
    """
    result = []             # will hold 1-letter list

    table = resDicStandard
    if nonstandard: table = resDic

## Python2.5    
##     table = resDic if nonstandard else resDicStandard
    if xtable:
        table = copy.copy( table )
        table.update( xtable )

    for aa in seq:
        try:
            aa = aa.lower()
            result +=  [ table[aa] ]
        except:
            result = result + [unknown]
    return result


def single2longAA( seq ):
    """
    Convert string of 1-letter AA code into list of 3-letter AA codes.
    
    :param seq: amino acid sequence in 1-letter code
    :type  seq: str
    
    :return: list with the amino acids in 3-letter code
    :rtype: [str]
    """
    ## invert AA dict
    invTab = {}

    for key in aaDicStandard:
        invTab[ aaDicStandard[key] ] = key

    result = []
    for aa in seq:
        try:
            aa = aa.upper()
            result += [ invTab[aa].upper() ]
        except:
            EHandler.warning("unknown residue: " + str(aa))
            result += ['Xaa']

    return result


def cmpAtoms( a1, a2 ):
    """
    Comparison function for bringing atoms into standard order
    within residues as defined by :class:`atomDic`.
    
    :param a1: atom dictionary
    :type  a1: CrossView or equivalent dictionary
    :param a2: atom dictionary
    :type  a2: CrossView or equivalent dictionary
    
    :return: int or list of matching positions
    :rtype: [-1|0|1]   
    """
    ## get standard order within residues
    target = atomDic[ a1['residue_name'] ]

    i1 = len( target )
    if a1['name'] in target:
        i1 = target.index( a1['name'] )

    i2 = len( target )
    if a2['name'] in target:
        i2 = target.index( a2['name'] )

    return (i1 > i2) - (i1 < i2)


def sortAtomsOfModel( model ):
    """
    Sort atoms within residues into the standard order defined in :class:`atomDic`.
    
    :param model: model to sort
    :type  model: PDBModel

    :return: model with sorted atoms
    :rtype: PDBModel
    """
    ## make a copy
    model = model.take( model.atomRange() )

    ## sort atoms
    model = model.sort( model.argsort( cmpAtoms ) )

    return model



#############
##  TESTING        
#############
from . import test as BT
        
class Test(BT.BiskitTest):
    """Test case"""

    def test_molUtils( self ):
        """molUtils test"""
        from biskit import PDBModel

        S = self
        
        ## load a structure
        S.m = PDBModel( t.testRoot('lig/1A19.pdb' ))
        S.model_1 = S.m.compress( S.m.maskProtein() )

        ## now sort in standard order
        S.model_2 = sortAtomsOfModel( S.model_1)

        ## compare the atom order
        cmp = []
        for a in S.model_1.atomRange():
            cmp += [ cmpAtoms( S.model_1.atoms[a], S.model_2.atoms[a] )]

        self.assertEqual( N0.sum(cmp), 159 )

        ## get the primaty sequence as a string
        S.seq = S.model_1.sequence()

        ## convert it to a list of three letter code
        S.seq=single2longAA(S.seq)

        ## convert it to a list in one letter code
        S.seq=singleAA(S.seq)

        self.assertEqual( ''.join(S.seq), S.model_1.sequence() )

if __name__ == '__main__':
    
    BT.localTest()
