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
Match 2 sequences against each other, deleting all positions that differ.
compareStructures() compares sequences of 2 structures and returns
a residue mask for each of them. 
"""

## see: https://www.python.org/dev/peps/pep-0366/
## allow relative imports when calling module as main script for testing
if __name__ == "__main__" and __package__ is None:
    import biskit
    __package__ = "biskit"

from .core import oldnumeric as N0
from . import tools as T
from .core.difflib_old import SequenceMatcher

import numpy as N


def getOpCodes( seq_1, seq_2 ):
    """
    Compares two sequences and returns a list with the information
    needed to convert the first one sequence into the second.

    :param seq_1: list of single letters
    :type seq_1: [ str ]
    :param seq_2: list of single letters
    :type seq_2: [ str ]
    
    :return: Optimization code from difflib::
             [('delete', 0, 1, 0, 0), ('equal', 1, 4, 0, 3),
              ('insert', 4, 4, 3, 4), ('equal', 4, 180, 4, 180)]
    :rtype: [tuples]
    """
    seqDiff = SequenceMatcher( None, ''.join(seq_1) , ''.join(seq_2) )
    seqDiff = seqDiff.get_opcodes()

    return seqDiff


def getSkipLists( seqDiff ):
    """
    Extracts information about what residues that have to be removed
    from sequence 1 (delete code) and sequence 2 (insert code).
    Returns deletion codes in the format (start_pos, length).

    :param seqDiff: opcodes
    :type  seqDiff: [tuples]

    :return: Lists of tuples containing regions of the sequences that
             should be deteted. Example::
               strucDel_1 = [(0, 1), (180, 4)]
               strucDel_2 = [(3, 1), (207, 4)]
             
    :rtype: [tuple], [tuple]
    """
    strucDel_2 = []
    strucDel_1 = []

    i=0
    for list in seqDiff:

        # residues to be deleted in sequence 1
        if str( seqDiff[i][0] ) == 'delete':
            strucDel_1 = strucDel_1 + [ (list[1],list[2]-list[1]) ]

        # residues to be deleted in sequence 2
        if str( seqDiff[i][0] ) == 'insert':
            strucDel_2 = strucDel_2 + [ (list[3],list[4]-list[3]) ]

        i += 1    

    return strucDel_1, strucDel_2



def getEqualLists( seqDiff ):
    """
    Extract information about regions in the sequences that are equal.
    Returns deletion codes in the format (start_pos, length).

    :param seqDiff: opcodes
    :type  seqDiff: [tuples]

    :return: Lists of tuples containing regions of the sequences that
             are equal. Example::
               strucEqual_1 = [(0, 216)]
               strucEqual_2 = [(0, 216)]
    :rtype: [tuple], [tuple] 
    """
    strucEqual_1 = []
    strucEqual_2 = []

    i=0
    for list in seqDiff:

        if str( seqDiff[i][0] ) == 'equal':
            strucEqual_1 = strucEqual_1 + [ (list[1],list[2]-list[1]) ]
            strucEqual_2 = strucEqual_2 + [ (list[3],list[4]-list[3]) ]

        i += 1

    return strucEqual_1, strucEqual_2



def expandRepeatsLeft( s, start, end, length=1 ):
    """recursively identify sequence repeats on left edge of s[start:end]"""
    core = s[start:end]

    if start-length>=0 and s[ start-length : start ] == core[0 : length]:
        start -= length
        start = expandRepeatsLeft( s, start, end )

    return start

def expandRepeatsRight( s, start, end, length=1 ):
    """recursively identify sequence repeats on right edge of s[start:end]"""
    core = s[start:end]

    if end+length<=len(s) and s[ end: end+length ] == core[-length:end]:
        end += length
        end = expandRepeatsRight( s, start, end, length )

    return end

def expandRepeats( s, start, size ):
    """
    Expand a text fragment within a larger string so that it includes any
    sequence repetitions to its right or left edge. 
    
    Example:
       ABC[BC]CCCDE  -> A[BCCCC]DE
       
    The idea here is to avoid alignment missmatches due to duplications. 
    The above to sequences could be aligned in several ways, for example:
       A--BC---DE      AB----C-DE
       ABCBCCCCDE  or  ABCBCCCCDE
    We don't know for sure which positions should be kept and which positions
    should be deleted in the longer string. So the most conservative approach
    is to remove the whole ambiguous fragment.
    
    :param s: input string
    :type  s: str
    :param start: start position of text fragment
    :type  start: int
    :param size: size of text fragment
    :type  size: int
    
    :return: start and size of expanded fragment
    :rtype: (int, int)
    """
    end = start + size
    left = [  expandRepeatsLeft(s,start,end,l) for l in range(size+1) ]
    right= [ expandRepeatsRight(s,start,end,l) for l in range(size+1) ]

    left = min(left)
    right= max(right)
    
    return left, right-left


def getEqual( seqAA, seqNr, equalList ):
    """
    Gets only the postions of the sequences that are equal according
    to the OpCodes. This should not be nessesary but might be usefull
    to skip 'replace' OpCode.

    :param seqAA: list with the amino acid sequence in one letter code
    :type  seqAA: [str]
    :param seqNr: list with the amino acid postitons
    :type  seqNr: [int]
    :param equalList: Lists of tuples containing regions of the sequences that
                      are equal
    :type  equalList: [tuple], [tuple]

  
    :return: lists of amino acids and positions where equal
    :rtype: [str], [int]  
    """
    equalSeqAA = []
    equalSeqNr = []

    # delete residues in delList
    for equal in equalList:

        equalSeqAA = equalSeqAA + seqAA[equal[0]:equal[0]+equal[1]] 
        equalSeqNr = equalSeqNr + seqNr[equal[0]:equal[0]+equal[1]] 

    return equalSeqAA, equalSeqNr


def del2mask( seq, *delpos ):
    """convert list of (from, to) delete positions into a mask of 0 or 1"""
    mask = N0.ones( len(seq) )

    for start, size in delpos:
        mask.put( range( start, start+size), 0 )
        
    return mask

def compareSequences( seqAA_1, seqAA_2 ):
    """
    """
    seqAA_1 = list( seqAA_1 )
    seqAA_2 = list( seqAA_2 )
    seqNr_1 = list(range( len( seqAA_1 )))  ## try removing list()
    seqNr_2 = list(range( len( seqAA_2 )))

    # get mask
    mask_1 = N0.zeros( len( seqNr_1 ) )
    mask_2 = N0.zeros( len( seqNr_2 ) )

    # compare sequences
    seqDiff = getOpCodes( seqAA_1, seqAA_2)

    # get delete lists
    del_1, del_2 =  getSkipLists( seqDiff )
    
    del_1 = [ expandRepeats( seqAA_1, *pos ) for pos in del_1 ]
    del_2 = [ expandRepeats( seqAA_2, *pos ) for pos in del_2 ]
    
    mask1 = del2mask( seqAA_1, *del_1 )
    mask2 = del2mask( seqAA_2, *del_2 )
 
    seqAA_1 = N0.compress( mask1, seqAA_1 ).tolist()
    seqNr_1 = N0.compress( mask1, seqNr_1 ).tolist()
    seqAA_2 = N0.compress( mask2, seqAA_2 ).tolist()
    seqNr_2 = N0.compress( mask2, seqNr_2 ).tolist()
    
    # get equal parts
    seqDiff = getOpCodes( seqAA_1, seqAA_2 )
    equal_1, equal_2 = getEqualLists( seqDiff ) 
    seqAA_1, seqNr_1 = getEqual( seqAA_1, seqNr_1, equal_1)
    seqAA_2, seqNr_2 = getEqual( seqAA_2, seqNr_2, equal_2 )

    N0.put( mask_1, seqNr_1 , 1 )
    N0.put( mask_2, seqNr_2 , 1 )

    return mask_1, mask_2

def compareModels( model_1, model_2 ):
    """
    Initiates comparison of the sequences of two structure objects and
    returns two equal sequence lists (new_seqAA_1 and new_seqAA_2 should
    be identical) and the corresponding residue position lists.

    :param model_1: model
    :type  model_1: PDBModel
    :param model_2: model
    :type  model_2: PDBModel

    :return: tuple of atom masks for model_1 and model_2::
              e.g. ( [0001011101111111], [1110000111110111] )
    :rtype: ([1|0...],[1|0...])
    """
    # get sequence AA and Nr strings
    seqAA_1 = model_1.sequence()
    seqAA_2 = model_2.sequence()

    return compareSequences( seqAA_1, seqAA_2 )


#############
##  TESTING        
#############
import biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test case"""

    def test_match2seq(self):
        """match2seq test"""
        ## Reading pdb files
        lig_traj = T.load( T.testRoot() + '/lig_pcr_00/traj.dat' )[:2]
        m = [ m.compress( m.maskProtein() ) for m in lig_traj ]

        ## make the models different
        m[1].removeRes(['ALA'])

        mask1, mask2 = compareModels( m[0], m[1] )

        if self.local:
            print('Reading and comparing two models')

            print('\nResidue masks to make the two maodels equal')
            print('mask1\n', mask1)
            print('mask2\n', mask2)            
            globals().update( locals() )

        self.assertTrue( N.all(mask1 == self.EXPECT[0] ) )
        self.assertTrue( N.all(mask2 == self.EXPECT[1] ) )
        
    def test_sequenceRepeats(self):
        """match2seq sequence repeat test"""
        seq1 = 'ABCDEFG~~~~~~~~~~~~~~~'
        seq2 = '~~~~~'
        mask1, mask2 = compareSequences( seq1, seq2 )
        self.assertTrue( N.all( mask1 == N0.zeros( len(seq1 ) )) )
        self.assertTrue( N.all( mask2 == N0.zeros( len(seq2 ) )) )


    EXPECT =  N0.array([1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                       1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                       1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1,
                       1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1,
                       1, 1, 1, 1],N0.Int),\
              N0.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],N0.Int)



if __name__ == '__main__':

    BT.localTest()
    



