# This module handles input and output of PDB files.
#
# Written by Konrad Hinsen <hinsen@cnrs-orleans.fr>
# Last revision: 2010-10-20
#
# reduced to bare-bones low-level PDB parsing class
# without dependency to numpy.oldnumeric by 
# Raik Gruenberg, 2016

"""
Parsing and writing of Protein Data Bank (PDB) files

This module provides classes that represent PDB (Protein Data Bank)
files and configurations contained in PDB files. It provides access to
PDB files on two levels: low-level (line by line) and high-level
(chains, residues, and atoms).

Caution: The PDB file format has been heavily abused, and it is
probably impossible to write code that can deal with all variants
correctly. This modules tries to read the widest possible range of PDB
files, but gives priority to a correct interpretation of the PDB
format as defined by the Brookhaven National Laboratory.

A special problem are atom names. The PDB file format specifies that
the first two letters contain the right-justified chemical element
name. A later modification allowed the initial space in hydrogen names
to be replaced by a digit. Many programs ignore all this and treat the
name as an arbitrary left-justified four-character name. This makes it
difficult to extract the chemical element accurately; most programs
write the '"CA"' for C_alpha in such a way that it actually stands for
a calcium atom. For this reason a special element field has been added
later, but only few files use it. In the absence of an element field,
the code in this module attempts to guess the element using all information
available.

The low-level routines in this module do not try to deal with the atom
name problem; they return and expect four-character atom names
including spaces in the correct positions. 

--Note (Raik)--

    The high-level classes have been removed from this version of the
    module. Refer to the original Scientific IO package

original documentation:

    The high-level routines use atom names without leading or trailing spaces,
    but provide and use the element field whenever possible. For output, they
    use the element field to place the atom name correctly, and for input, they
    construct the element field content from the atom name if no explicit
    element field is found in the file.

Except where indicated, numerical values use the same units and
conventions as specified in the PDB format description.


@undocumented: atom_format
@undocumented: anisou_format
@undocumented: conect_format
@undocumented: ter_format
@undocumented: model_format
@undocumented: header_format
@undocumented: cryst1_format
@undocumented: scalen_format
@undocumented: mtrixn_format
@undocumented: generic_format
@undocumented: export_filters
@undocumented: DummyChain
"""

## see: https://www.python.org/dev/peps/pep-0366/
## allow relative imports when calling module as main script for testing
if __name__ == "__main__" and __package__ is None:
    import biskit
    __package__ = "biskit.core.scientificIO"

from .TextFile import TextFile
from .FortranFormat import FortranFormat, FortranLine
from .PDBExportFilters import export_filters

import numpy as N
import copy, string

#
# Fortran formats for PDB entries
#
atom_format = FortranFormat('A6,I5,1X,A4,A1,A4,A1,I4,A1,3X,3F8.3,2F6.2,' +
                            '6X,A4,2A2')
anisou_format = FortranFormat('A6,I5,1X,A4,A1,A4,A1,I4,A1,1X,6I7,2X,A4,2A2')
conect_format = FortranFormat('A6,11I5')
ter_format = FortranFormat('A6,I5,6X,A4,A1,I4,A1')
model_format = FortranFormat('A6,4X,I4')
header_format = FortranFormat('A6,4X,A40,A9,3X,A4')
cryst1_format = FortranFormat('A6,3F9.3,3F7.2,1X,A11,I4')
scalen_format = FortranFormat('A6,4X,3F10.6,5X,F10.5')
mtrixn_format = FortranFormat('A6,1X,I3,3F10.6,5X,F10.5,4X,I1')
generic_format = FortranFormat('A6,A74')

#
# Amino acid and nucleic acid residues
#
amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'CYX', 'GLN', 'GLU', 'GLY',
               'HIS', 'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP', 'ILE', 'LEU',
               'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
               'ACE', 'NME', 'NHE']

nucleic_acids = [ 'A',  'C',  'G',  'I',  'T',  'U',
                 '+A', '+C', '+G', '+I', '+T', '+U',
                  'RA',  'RC',  'RG',  'RU',
                  'DA',  'DC',  'DG',  'DT',
                  'RA5',  'RC5',  'RG5',  'RU5',
                  'DA5',  'DC5',  'DG5',  'DT5',
                  'RA3',  'RC3',  'RG3',  'RU3',
                  'DA3',  'DC3',  'DG3',  'DT3',
                  'RAN',  'RCN',  'RGN',  'RUN',
                  'DAN',  'DCN',  'DGN',  'DTN',
                  ]

def defineAminoAcidResidue(symbol):
    """
    Make the parser recognize a particular residue type as an amino
    acid residue
    
    :param symbol: the three-letter code for an amino acid
    :type symbol: ``str``
    """
    symbol = symbol.upper()
    if symbol not in amino_acids:
        amino_acids.append(symbol)

def defineNucleicAcidResidue(symbol):
    """
    Make the parser recognize a particular residue type as an nucleic
    acid residue
    
    :param symbol: the one-letter code for a nucleic acid
    :type symbol: ``str``
    """
    symbol = symbol.upper()
    if symbol not in nucleic_acids:
        nucleic_acids.append(symbol)


#
# Low-level file object. It represents line contents as Python dictionaries.
# For output, there are additional methods that generate sequence numbers
# for everything.
#
class PDBFile:

    """
    X{PDB} file with access at the record level

    The low-level file access is handled by the module
    :class:`ScientificIO.TextFile`, therefore compressed files and URLs
    (for reading) can be used as well.
    """

    def __init__(self, file_or_filename, mode = 'r', subformat = None):
        """
        :param file_or_filename: the name of the PDB file, or a file object
        :type file_or_filename: ``str`` or C{file}
        :param mode: the file access mode, 'r' (read) or 'w' (write)
        :type mode: ``str``
        :param subformat: indicates a specific dialect of the PDB format.
                          Subformats are defined in
                          :class:`ScientificIO.PDBExportFilters`; they are used
                          only when writing.
        :type subformat: ``str`` or C{NoneType}
        """
        if isinstance(file_or_filename, str):
            self.file = TextFile(file_or_filename, mode)
        else:
            self.file = file_or_filename
        self.output = str.lower(mode[0]) == 'w'
        self.export_filter = None
        if subformat is not None:
            export = export_filters.get(subformat, None)
            if export is not None:
                self.export_filter = export()
        self.open = 1
        if self.output:
            self.data = {'serial_number': 0,
                         'residue_number': 0,
                         'chain_id': '',
                         'segment_id': ''}
            self.het_flag = 0
            self.chain_number = -1

    def readLine(self):
        """
        Return the contents of the next non-blank line (= record) The
        return value is a tuple whose first element (a string)
        contains the record type. For supported record types (HEADER,
        CRYST1, SCALEn, MTRIXn, ATOM, HETATM, ANISOU, TERM, MODEL,
        CONECT), the items from the remaining fields are put into a
        dictionary which is returned as the second tuple element. Most
        dictionary elements are strings or numbers; atom positions are
        returned as a vector, and anisotropic temperature factors are
        returned as a rank-2 tensor, already multiplied by 1.e-4.
        White space is stripped from all strings except for atom
        names, whose correct interpretation can depend on an initial
        space. For unsupported record types, the second tuple element
        is a string containing the remaining part of the record.

        :returns: the contents of one PDB record
        :rtype: ``tuple``
        """
        while 1:
            line = self.file.readline()
            if not line: return ('END','')
            if line[-1] == '\n': line = line[:-1]
            line = str.strip(line)
            if line: break
        line = str.ljust(line, 80)
        type = str.strip(line[:6])
        if type == 'ATOM' or type == 'HETATM':
            line = FortranLine(line, atom_format)
            data = {'serial_number': line[1],
                    'name': line[2],
                    'alternate': str.strip(line[3]),
                    'residue_name': str.strip(line[4]),
                    'chain_id': str.strip(line[5]),
                    'residue_number': line[6],
                    'insertion_code': str.strip(line[7]),
##                    'position': N.array(line[8:11]),
                    'position': line[8:11],
                    'occupancy': line[11],
                    'temperature_factor': line[12],
                    'segment_id': str.strip(line[13]),
                    'element': str.strip(line[14]),
                    'charge': str.strip(line[15])}
            return type, data
        elif type == 'ANISOU':
            line = FortranLine(line, anisou_format)
            data = {'serial_number': line[1],
                    'name': line[2],
                    'alternate': str.strip(line[3]),
                    'residue_name': str.strip(line[4]),
                    'chain_id': str.strip(line[5]),
                    'residue_number': line[6],
                    'insertion_code': str.strip(line[7]),
                    'u': 1.e-4* N.array( [[line[8], line[11], line[12]],
                                          [line[11], line[9] , line[13]],
                                          [line[12], line[13], line[10]]] ),
                    'segment_id': str.strip(line[14]),
                    'element': str.strip(line[15]),
                    'charge': str.strip(line[16])}
            return type, data
        elif type == 'TER':
            line = FortranLine(line, ter_format)
            data = {'serial_number': line[1],
                    'residue_name': str.strip(line[2]),
                    'chain_id': str.strip(line[3]),
                    'residue_number': line[4],
                    'insertion_code': str.strip(line[5])}
            return type, data
        elif type == 'CONECT':
            line = FortranLine(line, conect_format)
            data = {'serial_number': line[1],
                    'bonded': [i for i in line[2:6] if i > 0],
                    'hydrogen_bonded': [i for i in line[6:10] if i > 0],
                    'salt_bridged': [i for i in line[10:12] if i > 0]}
            return type, data
        elif type == 'MODEL':
            line = FortranLine(line, model_format)
            data = {'serial_number': line[1]}
            return type, data
        elif type == 'HEADER':
            line = FortranLine(line, header_format)
            data = {'compound': line[1],
                    'date': line[2],
                    'pdb_code': line[3]}
            return type, data
        elif type == 'CRYST1':
            line = FortranLine(line, cryst1_format)
            data = {'a': line[1],
                    'b': line[2],
                    'c': line[3],
                    'alpha': line[4],
                    'beta': line[5],
                    'gamma': line[6],
                    'space_group': line[7],
                    'z': line[8]}
            return type, data
        elif type[:-1] == 'SCALE':
            line = FortranLine(line, scalen_format)
            data = {'s1': line[1],
                    's2': line[2],
                    's3': line[3],
                    'u': line[4]}
            return type, data
        elif type[:-1] == 'MTRIX':
            line = FortranLine(line, mtrixn_format)
            data = {'serial': line[1],
                    'm1': line[2],
                    'm2': line[3],
                    'm3': line[4],
                    'v': line[5],
                    'given': line[6] == 1}
            return type, data
        else:
            return type, line[6:]

    def writeLine(self, type, data):
        """
        Write a line using record type and data dictionary in the
        same format as returned by readLine(). Default values are
        provided for non-essential information, so the data dictionary
        need not contain all entries.

        :param type: PDB record type
        :type type: ``str``
        :param data: PDB record data
        :type data: ``tuple``
        """
        if self.export_filter is not None:
            type, data = self.export_filter.processLine(type, data)
            if type is None:
                return
        line = [type]
        if type == 'ATOM' or type == 'HETATM':
            format = atom_format
            position = data['position']
            line = line + [data.get('serial_number', 1),
                           data.get('name'),
                           data.get('alternate', ''),
                           str.rjust(data.get('residue_name', ''), 3),
                           data.get('chain_id', ''),
                           data.get('residue_number', 1),
                           data.get('insertion_code', ''),
                           position[0], position[1], position[2],
                           data.get('occupancy', 0.),
                           data.get('temperature_factor', 0.),
                           data.get('segment_id', ''),
                           str.rjust(data.get('element', ''), 2),
                           data.get('charge', '')]
        elif type == 'ANISOU':
            format = anisou_format
            u = 1.e4*data['u']
            u = [int(u[0,0]), int(u[1,1]), int(u[2,2]),
                 int(u[0,1]), int(u[0,2]), int(u[1,2])]
            line = line + [data.get('serial_number', 1),
                           data.get('name'),
                           data.get('alternate', ''),
                           str.rjust(data.get('residue_name'), 3),
                           data.get('chain_id', ''),
                           data.get('residue_number', 1),
                           data.get('insertion_code', '')] \
                        + u \
                        + [data.get('segment_id', ''),
                           str.rjust(data.get('element', ''), 2),
                           data.get('charge', '')]
        elif type == 'TER':
            format = ter_format
            line = line + [data.get('serial_number', 1),
                           str.rjust(data.get('residue_name'), 3),
                           data.get('chain_id', ''),
                           data.get('residue_number', 1),
                           data.get('insertion_code', '')]
        elif type == 'CONECT':
            format = conect_format
            line = line + [data.get('serial_number')]
            line = line + (data.get('bonded', [])+4*[None])[:4]
            line = line + (data.get('hydrogen_bonded', [])+4*[None])[:4]
            line = line + (data.get('salt_bridged', [])+2*[None])[:2]
        elif type == 'MODEL':
            format = model_format
            line = line + [data.get('serial_number')]
        elif type == 'CRYST1':
            format = cryst1_format
            line = line + [data.get('a'), data.get('b'), data.get('c'),
                           data.get('alpha'), data.get('beta'),
                           data.get('gamma'),
                           data.get('space_group'),
                           data.get('z')]
        elif type[:-1] == 'SCALE':
            format = scalen_format
            line = line + [data.get('s1'), data.get('s2'), data.get('s3'),
                           data.get('u')]
        elif type[:-1] == 'MTRIX':
            format = scalen_format
            line = line + [data.get('serial'),
                           data.get('m1'), data.get('m2'), data.get('m3'),
                           data.get('v'), int(data.get('given'))]
        elif type == 'HEADER':
            format = header_format
            line = line + [data.get('compound', ''), data.get('date', ''),
                           data.get('pdb_code')]
        else:
            format = generic_format
            line = line + [data]
        self.file.write(str(FortranLine(line, format)) + '\n')

    def writeComment(self, text):
        """
        Write text into one or several comment lines.
        Each line of the text is prefixed with 'REMARK' and written
        to the file.

        :param text: the comment contents
        :type text: ``str``
        """
        while text:
            eol = str.find(text,'\n')
            if eol == -1:
                eol = len(text)
            self.file.write('REMARK %s \n' % text[:eol])
            text = text[eol+1:]

    def writeAtom(self, name, position, occupancy=0.0, temperature_factor=0.0,
                  element='', alternate=''):
        """
        Write an ATOM or HETATM record using the information supplied.
        The residue and chain information is taken from the last calls to
        the methods :class:`nextResidue` and :class:`nextChain`.

        :param name: the atom name
        :type name: ``str``
        :param position: the atom position
        :type position: :class:`numpy.ndarray`
        :param occupancy: the occupancy
        :type occupancy: ``float``
        :param temperature_factor: the temperature factor (B-factor)
        :type temperature_factor: ``float``
        :param element: the chemical element
        :type element: ``str``
        :param alternate: the alternate location character
        :type element: ``str``
        """
        if self.het_flag:
            type = 'HETATM'
        else:
            type = 'ATOM'
        name = str.upper(name)
        if element != '' and len(element) == 1 and name and name[0] == element and len(name) < 4:
            name = ' ' + name
        self.data['name'] = name
        self.data['position'] = position
        self.data['serial_number'] = (self.data['serial_number'] + 1) % 100000
        self.data['alternate'] = alternate
        self.data['occupancy'] = occupancy
        self.data['temperature_factor'] = temperature_factor
        self.data['element'] = element
        self.writeLine(type, self.data)

    def nextResidue(self, name, number = None, terminus = None):
        """
        Signal the beginning of a new residue, starting with the
        next call to :class:`writeAtom`.

        :param name: the residue name
        :type name: ``str``
        :param number: the residue number. If ``None``, the residues
                       will be numbered sequentially, starting from 1.
        :type number: ``int`` or C{NoneType}
        :param terminus: ``None``, "C", or "N". This information
                         is passed to export filters that can use this
                         information in order to use different atom or
                         residue names in terminal residues.
        """
        name  = str.upper(name)
        if self.export_filter is not None:
            name, number = self.export_filter.processResidue(name, number,
                                                             terminus)
        self.het_flag =  not (name in amino_acids or name in nucleic_acids)
        self.data['residue_name'] = name
        self.data['residue_number'] = (self.data['residue_number'] + 1) % 10000
        self.data['insertion_code'] = ''
        if number is not None:
            if isinstance(number, int):
                if number >= 0:
                    self.data['residue_number'] = number % 10000
                else:
                    self.data['residue_number'] = -((-number) % 1000)
            else:
                self.data['residue_number'] = number.number % 10000
                self.data['insertion_code'] = number.insertion_code

    def nextChain(self, chain_id = None, segment_id = ''):
        """
        Signal the beginning of a new chain.

        :param chain_id: a chain identifier. If ``None``, consecutive letters
                         from the alphabet are used.
        :type chain_id: ``str`` or C{NoneType}
        :param segment_id: a chain identifier
        :type segment_id: ``str``
        """
        if chain_id is None:
            self.chain_number = (self.chain_number + 1) % len(self._chain_ids)
            chain_id = self._chain_ids[self.chain_number]
        if self.export_filter is not None:
            chain_id, segment_id = \
                      self.export_filter.processChain(chain_id, segment_id)
        self.data['chain_id'] = (chain_id+' ')[:1]
        self.data['segment_id'] = (segment_id+'    ')[:4]
        self.data['residue_number'] = 0

    _chain_ids = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    def terminateChain(self):
        """
        Signal the end of a chain.
        """
        if self.export_filter is not None:
            self.export_filter.terminateChain()
        self.data['serial_number'] = (self.data['serial_number'] + 1) % 100000
        self.writeLine('TER', self.data)
        self.data['chain_id'] = ''
        self.data['segment_id'] = ''
        
    def close(self):
        """
        Close the file. This method **must** be called for write mode
        because otherwise the file will be incomplete.
        """
        if self.open:
            if self.output:
                self.file.write('END\n')
            self.file.close()
            self.open = 0

    def __del__(self):
        try:
            self.close()
        except:
            pass

############
## TESTING
############

import biskit.test as BT

class Test(BT.BiskitTest):
    """Test case"""

    def prepare(self):
        import tempfile
        self.f_pdbcopy = tempfile.mktemp(suffix='_test_PDB_copy.pdb')
        
    def cleanUp(self):
        import biskit.tools as T
        T.tryRemove(self.f_pdbcopy)

    def test_PDB( self ):
        """PDB with Anisotropy read and write"""
        import biskit.tools as T
        import os.path as osp        
        
        fname = osp.join(T.testRoot(subfolder='rec'),'1A2P_rec_original.pdb')
        f = PDBFile( fname )
        fcopy = PDBFile(self.f_pdbcopy, mode='w', subformat='xplor')

        while 1:
            type, data = f.readLine()
            fcopy.writeLine(type, data)
            if type == 'END':
                break
        fcopy.close()

if __name__ == '__main__':
    
    BT.localTest()
