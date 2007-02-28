##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2006 Raik Gruenberg & Johan Leckner
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
## $Revision$
## last $Date$
## last $Author$
"""
Match 2 sequences against each other, deleting all positions that differ.
compareStructures() compares sequences of 2 structures and returns
a residue mask for each of them. 
"""

import Numeric as N
import Biskit.tools as T
import sys
from Biskit.difflib_old import SequenceMatcher
#from difflib import SequenceMatcher


def sequenceList( model ):
    """
    Extracts a one letter amino acid sequence list
    
    @param model: model
    @type  model: PDBModel
    
    @return: sequence ['A','G','R','T',.....]
    @rtype: [str]
    """
    return [ i for i in model.sequence() ]


def getOpCodes( seq_1, seq_2 ):
    """
    Compares two sequences and returns a list with the information
    needed to convert the first one sequence into the second.

    @param seq_1: list of single letters
    @type seq_1: [ str ]
    @param seq_2: list of single letters
    @type seq_2: [ str ]
    
    @return: Optimization code from difflib::
             [('delete', 0, 1, 0, 0), ('equal', 1, 4, 0, 3),
              ('insert', 4, 4, 3, 4), ('equal', 4, 180, 4, 180)]
    @rtype: [tuples]
    """
    seqDiff = SequenceMatcher( None, ''.join(seq_1) , ''.join(seq_2) )
    seqDiff = seqDiff.get_opcodes()

    return seqDiff


def getSkipLists( seqDiff ):
    """
    Extracts information about what residues that have to be removed
    from sequence 1 (delete code) and sequence 2 (insert code).
    Returns deletion codes in the format (start_pos, length).

    @param seqDiff: opcodes
    @type  seqDiff: [tuples]

    @return: Lists of tuples containing regions of the sequences that
             should be deteted. Example::
               strucDel_1 = [(0, 1), (180, 4)]
               strucDel_2 = [(3, 1), (207, 4)]
             
    @rtype: [tuple], [tuple]
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

    @param seqDiff: opcodes
    @type  seqDiff: [tuples]

    @return: Lists of tuples containing regions of the sequences that
             are equal. Example::
               strucEqual_1 = [(0, 216)]
               strucEqual_2 = [(0, 216)]
    @rtype: [tuple], [tuple] 
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



def repeateInMatch( deletedRes ):
    """
    Takes a string and returns a list of strings containing internal repeates.

    @param deletedRes: string of one letter amino acid sequence
    @type  deletedRes: str
    
    @return: smallest self repete. Example::
                'GSGS'     returns  'GS'
                'ABABABAB'    "     'ABAB' and 'AB'
                'FFF'         "     'F'
    @rtype: str
    """
    segment = []
    repete = []
    word = []
    checkList = []

    # the maximun length of a word is half of the deleted residues
    length = len( deletedRes )
    half = range( 1, length/2+1 )

    if len( half ) > 1:
        half.reverse()

    # get a list of the number of possible repetes
    #
    # example: for a string 'ABABABAB' (lengt 8) return [2,4,8] corresponding
    #          to 'ABAB','ABAB' or 'AB','AB',AB','AB' or 'A','B','A',...
    for i in half:

        if length/(i*1.0) - int( length/(i*1.0) ) == 0.0:
            segment.append(length/i)

    # get the corresponding words
    # example: 'ABABABAB' and seg = 1 -> ['ABAB','ABAB']
    for seg in segment:

        wordLen = length/seg
        i=0
        word = []

        while i < length:
            word = word + [deletedRes[i:i+wordLen]]
            i = i + wordLen

        # compare the identity of all words of length i
        # to all other words of the same length
        j=1
        k=2
        score = 0

        while j <= seg:

            while k <= seg:

                if word[j-1] == word[k-1]:
                    score += 1
                k += 1

            k = j +2
            j += 1
        # if all words are the same, add to list of strings to return
        if score == (seg*(seg-1))/2:

            checkList = checkList + [word[0]]

    return checkList    



def delete( seqAA, seqNr, delList ):
    """
    1. Takes a amino acid and a sequence list and deletes positions
       according to the information given in the delList.

    2. Furthermore, compares the deleted sequence with the following
       and preceding sequence. If they are identical also these
       residues are deleted.

    3. If the sequence contains internal repeates (as in the sequence
       'GSGS') and 2) does not apply the preceding and following
       sequence is also scanned for this sequence (here 'GS').

    @param seqAA: list with the amino acid sequence in one letter code
    @type  seqAA: [str]
    @param seqNr: list with the amino acid postitons
    @type  seqNr: [int]
    @param delList: list of residues to be deleted (postiton, length)
    @type  delList: [tuple]

    @return: sequence and positions, example::
               seqAA - ['A','G','R','T',.....]
               seqNr - [ 0,  1,  2,  3 ,.....]  
               delList - [(0, 2), (180, 4)]
               ->  seqAA - ['R','T',.....]
                   seqNr - [ 2,  3 ,.....]
                   
    @rtype: [str],[int]
    """

    offset = 0 # keep track of the number of deteted residues

    # iterate over the lists in delList
    for list in delList:

        deletedRes = ''     # seq of deleted residues ('GR')
        followingBlock = '' # seq of following res 
        precedingBlock = '' # seq of preceding res

        start = list[0]-offset 
        length = list[1]
        positions = range( start, start+length )

        # collect deteted residue sequence in a string
        # an the preceding and following sequence of same length
        for pos in positions:

            deletedRes = deletedRes + seqAA[ pos ]
            try:
                followingBlock = followingBlock  + seqAA[ pos + length ]
            except:
                pass
##                sys.stderr.write("N'del ")

            try:
                precedingBlock = precedingBlock + seqAA[ pos - length ]
            except:
                pass
##                sys.stderr.write("C'del ")

        # get a list of internal delete sequence repeats
        intraRepeteList = repeateInMatch(deletedRes)

        to_be_removed = [] 
        for pos in positions:

            to_be_removed.append( pos )
            offset += 1

            # if the sequence block of delList matches the following block
            # delete it as well
            if deletedRes == followingBlock:

                to_be_removed.append( pos+length )
                offset += 1

            # ... and the same goes for the previous block
            if deletedRes == precedingBlock:

                to_be_removed.append( pos-length )
                offset += 1

        # if the preceding block was not deleted, check for preceding
        # intra repetes
        hit = 0  # to be able to stop if hit is found
        if deletedRes != precedingBlock:

            while hit == 0:

                for repete in intraRepeteList:

                    repeteSeq = N.sum( seqAA[ positions[0]-len(repete) : positions[0] ] )
                    repetePos = range( positions[0]-len(repete), positions[0] )

                    if repete == repeteSeq:

                        hit += 1

                        for pos in repetePos:

                            to_be_removed.append( pos )
                            offset += 1
                hit += 1

        # ... and following intra repetes
        hit = 0
        if deletedRes != followingBlock:

            while hit == 0:

                for repete in intraRepeteList:

                    repeteSeq = N.sum( seqAA[ positions[0]+length : positions[0]+length+len(repete) ] )
                    repetePos =  range( positions[0]-length, positions[0]+length+len(repete) )

                    if repete == repeteSeq:

                        hit += 1

                        for pos in repetePos:

                            to_be_removed.append( pos )
                        offset += 1
                hit += 1

        # sort and inverse list of resides to be removed
        if len(to_be_removed) > 1:

            to_be_removed.sort()
            to_be_removed.reverse()

        # delete from the end and forward
        for r in to_be_removed:

            del seqAA[r]
            del seqNr[r]


    return seqAA, seqNr



def getEqual( seqAA, seqNr, equalList ):
    """
    Gets only the postions of the sequences that are equal according
    to the OpCodes. This should not be nessesary but might be usefull
    to skip 'replace' OpCode.

    @param seqAA: list with the amino acid sequence in one letter code
    @type  seqAA: [str]
    @param seqNr: list with the amino acid postitons
    @type  seqNr: [int]
    @param equalList: Lists of tuples containing regions of the sequences that
                      are equal
    @type  equalList: [tuple], [tuple]

  
    @return: lists of amino acids and positions where equal
    @rtype: [str], [int]  
    """
    equalSeqAA = []
    equalSeqNr = []

    # delete residues in delList
    for equal in equalList:

        equalSeqAA = equalSeqAA + seqAA[equal[0]:equal[0]+equal[1]] 
        equalSeqNr = equalSeqNr + seqNr[equal[0]:equal[0]+equal[1]] 


    return equalSeqAA, equalSeqNr


def iterate( seqAA_1, seqNr_1, seqAA_2, seqNr_2, del_1, del_2 ):
    """
    Delete residues until no more deletions are indicated by the
    sequence matcher. Return the final sequences

    @param seqAA_1: list with the amino acid sequence in one letter code
    @type  seqAA_1: [str]    
    @param seqNr_1: list with the amino acid postitons
    @type  seqNr_1: [int]
    @param seqAA_2: list with the amino acid sequence in one letter code
    @type  seqAA_2: [str]    
    @param seqNr_2: list with the amino acid postitons
    @type  seqNr_2: [int]
    @param del_1: Lists of tuples containing regions of the sequences that
                  should be deteted
    @type  del_1: [tuple]
    @param del_2: Lists of tuples containing regions of the sequences that
                  should be deteted
    @type  del_2: [tuple]
    
    @return: the final sequence and position lists
    @rtype: [str], [int], [str], [int]
    """
    while del_1 <> [] or del_2 <> []:
        seqAA_1, seqNr_1  = delete(seqAA_1, seqNr_1, del_1)
        seqAA_2, seqNr_2  = delete(seqAA_2, seqNr_2, del_2)

        seqDiff = getOpCodes( seqAA_1, seqAA_2 )
        del_1, del_2 = getSkipLists( seqDiff )

    return  seqAA_1, seqNr_1, seqAA_2, seqNr_2


def compareModels( model_1, model_2 ):
    """
    Initiates comparison of the sequences of two structure objects and
    returns two equal sequence lists (new_seqAA_1 and new_seqAA_2 should
    be identical) and the corresponding residue position lists.

    @param model_1: model
    @type  model_1: PDBModel
    @param model_2: model
    @type  model_2: PDBModel

    @return: tuple of atom masks for model_1 and model_2::
              e.g. ( [0001011101111111], [1110000111110111] )
    @rtype: ([1|0...],[1|0...])
    """
    # get sequence AA and Nr strings
    seqAA_1 = sequenceList( model_1 )
    seqAA_2 = sequenceList( model_2 )
    seqNr_1 = range( len( seqAA_1 ) )
    seqNr_2 = range( len( seqAA_2 ) )

    # get mask
    mask_1 = N.zeros( len( seqNr_1 ), N.Int )
    mask_2 = N.zeros( len( seqNr_2 ), N.Int )

    # compare sequences
    seqDiff = getOpCodes( seqAA_1, seqAA_2)

    # get delete lists
    del_1, del_2 =  getSkipLists( seqDiff )

    # delete residues until del_1 = del_2 = []
    seqAA_1, seqNr_1, seqAA_2, seqNr_2 = \
                 iterate( seqAA_1, seqNr_1, seqAA_2, seqNr_2, del_1, del_2 )

    # get equal parts
    seqDiff = getOpCodes( seqAA_1, seqAA_2 )
    equal_1, equal_2 = getEqualLists( seqDiff ) 
    seqAA_1, seqNr_1 = getEqual( seqAA_1, seqNr_1, equal_1)
    seqAA_2, seqNr_2 = getEqual( seqAA_2, seqNr_2, equal_2 )

    N.put( mask_1, seqNr_1 , 1 )
    N.put( mask_2, seqNr_2 , 1 )

    return mask_1, mask_2



#############
##  TESTING        
#############
import Biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test case"""

    def test_match2seq(self):
	"""match2seq test"""
        ## Reading pdb files
        lig_traj = T.Load( T.testRoot() + '/lig_pcr_00/traj.dat' )[:2]
        m = [ m.compress( m.maskProtein() ) for m in lig_traj ]

        ## make the models different
        m[1].removeRes(['ALA'])

        mask1, mask2 = compareModels( m[0], m[1] )

        if self.local:
            print 'Reading and comparing two models'

            print '\nResidue masks to make the two maodels equal'
            print 'mask1\n', mask1
            print 'mask2\n', mask2            
            globals().update( locals() )

        self.assertEqual( (mask1, mask2), self.EXPECT )


    EXPECT =  N.array([1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		       1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		       1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1,
		       1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1,
		       1, 1, 1, 1],N.Int),\
	      N.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],N.Int)



if __name__ == '__main__':

    BT.localTest()
    



