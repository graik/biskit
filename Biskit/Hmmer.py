## Class conservation
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
## $Revision$
## last $Author$
## last $Date$

"""
Search Hmmer Pfam database and retrieve conservation data.
"""

import tempfile
import os, re, string
import types
import mathUtils as math
import Numeric as N
import molUtils
import settings
import Biskit.Mod.modUtils as MU

## executables
hmmpfamExe = settings.hmmpfam_bin
hmmfetchExe = settings.hmmfetch_bin
hmmalignExe = settings.hmmalign_bin
hmmindexExe = settings.hmmindex_bin
hmmDatabase = settings.hmm_db


def hmmEmm2Prob( nullEmm, emmScore ):
    """
    Convert HMM profile emmisiion scores into emmission probabilities

    @param nullEmm: null scores
    @type  nullEmm: array
    @param emmScore: emmission scores
    @type  emmScore: array

    @return: null and emmission probabilities, for each amino acid
             in each position
    @rtype:  array( len_seq x 20 ), array( 1 x 20 )    
    """
    ## Null probabilities: prob = 2 ^ (nullEmm / 1000) * 1/len(alphabet)
    nullProb = N.power( 2, N.array( nullEmm )/1000.0 )*(1./20)

    ## Emmission probabilities: prob = nullProb 2 ^ (nullEmm / 1000)
    emmProb = nullProb * N.power( 2, ( emmScore/1000.0) )

    return emmProb, nullProb


def entropy( emmProb, nullProb ):
    """
    calculate entropy for normalized probabilities scaled by aa freq.
    emmProb & nullProb is shape 1,len(alphabet)

    @param emmProb: emmission probabilities
    @type  emmProb: array
    @param nullProb: null probabilities
    @type  nullProb: array

    @return: entropy value
    @rtype:  float
    """
    ## remove zeros to avoid log error
    emmProb = N.clip(emmProb, 1.e-10, 1.)

    return N.sum( emmProb * N.log(emmProb/nullProb) )


class HmmerError( Exception ):
    pass


class Hmmer:
    """
    Search Hmmer Pfam database and retrieve conservation score for model
    """

    def __init__(self, hmmdb = hmmDatabase ): 
        """
        @param hmmdb: Pfam hmm database
        @type  hmmdb: str
        """
        self.hmmdb = hmmdb
        self.tempDir = settings.tempDirShared

        self.fastaID = ''

        tempfile.tempdir = self.tempDir
        self.hmmFile = tempfile.mktemp('.hmm')
        self.fastaFile = tempfile.mktemp('.fasta')
        self.sub_fastaFile = tempfile.mktemp('_sub.fasta')


    def checkHmmdbIndex( self ):
        """
        Checks if the hmm database has been indexed or not,
        if not indexing will be done.
        """
        if not os.path.exists(self.hmmdb+'.ssi'):
            print 'HMMINDEX: Indexing hmm database. This will take a while'
            os.system(hmmindexExe + ' ' + self.hmmdb)


    def fasta( self, m , start=0, stop=None ):
        """
        Extract fasta sequence from model.

        @param m: model
        @type  m: PDBModel
        @param start: first residue
        @type  start: int
        @param stop: last residue
        @type  stop: int
        
        @return: fasta formated sequence and PDB code
        @rtype: string
        """
        if not stop:
            stop = m.lenResidues()

        s = m.sequence()[start:stop]

        n_chunks = len( s ) / 80

        result = ">%s \n" % m.pdbCode

        for i in range(0, n_chunks+1):

            if i * 80 + 80 < len( s ):
                chunk = s[i * 80 : i * 80 + 80]
            else:
                chunk = s[i * 80 :]

            result += chunk
        return  result, m.pdbCode


    def searchHmmdb( self, target, noSearch=None ):
        """
        Search hmm database with a sequence in fasta format.
        If the profile names have been provided - skip the search and
        only write the temporary sequence files.
          
        @param target: sequence 
        @type  target: PDBModel or fasta file
        @param noSearch: don't perform a seach
        @type  noSearch: 1 OR None
         
        @return: dictionary witn profile names as keys and a list of
                 lists containing information about the range where the
                 profile matches the sequence
        @rtype: dict, [list]
        """
        ## if target is a PDBModel
        if type(target) == types.InstanceType:
            fastaSeq, self.fastaID = self.fasta( target )

            ## write fasta sequence file
            fName = self.fastaFile
            seq = open( fName, 'w' )
            seq.write( fastaSeq )
            seq.close()

        ## else assume it is a fasta sequence file   
        else:
            if MU.verify_fasta(target):
                fName = target

        if noSearch:
            print 'Profiles provided - No search will be performed.'
            return None

        else:
            ## find matching hmm profiles
            print '\nSearching hmm database, this will take a while...'
            out = os.popen( hmmpfamExe + ' ' + self.hmmdb + ' ' + fName )

            matches = {}
            hits = []
            try:
                while 1:
                    l = out.readline()
                    ## get names and descriptions of matching profiles
                    if re.match('^-{8}\s{7,8}-{11}.+-{3}$', l):
                        m = string.split( out.readline() )
                        while len(m) != 0:
                            matches[m[0]] =  m[1:] 
                            m = string.split( out.readline() )

                    ## get hits, scores and alignment positions
                    if re.match('^-{8}\s{7,8}-{7}\s-{5}\s-{5}.+-{7}$', l):
                        h = string.split( out.readline() )
                        while len(h) != 0:
                            hits += [ h ] 
                            h = string.split( out.readline() )
                        break

            except:
                print '\nERROR parsing hmmpfam search result.'

            return matches, hits


    def selectMatches( self, matches, hits, score_cutoff=60 ,
                       eValue_cutoff = 1e-8 ):
        """
        select what hmm profiles to use based on score and e-Value cutoff

        @param matches: output from L{searchHmmdb}
        @type  matches: dict
        @param hits:  output from L{searchHmmdb}
        @type  hits: [list] 
        @param score_cutoff: cutoff value for an acceptable score
        @type  score_cutoff: float
        @param eValue_cutoff: cutoff value for an acceptable e-value
        @type  eValue_cutoff: float

        @return: dictionary with good matches {hmm_name : [[start,stop],[..]]}
        @rtype: dict     
        """
        ## Print list of matching profiles
        if len(matches) == 0:
            print '\nNo matching profiles found'
            return None
        if len(matches) == 1:
            print '\nFound only one matching profile.'
        if len(matches) > 1:
            print '\nFound more than one hit.\n'

        try:
            print '\tName -- Nr hits -- Score -- E-value -- Description'
            for k in matches.keys():
                print '\t', k, ' --', matches[k][-1], ' -- ',\
                      matches[k][-3], ' -- ', matches[k][-2],' -- ',\
                      string.join(matches[k][0:-3])
        except:
            pass

        ## select to use profiles from searchHits  given the given criteria
        hmmHits = {}

        for h in range( len(hits) ):

            hmmName = hits[h][0]
            score = float( hits[h][-2] )
            eValue = float( hits[h][-1] )

            ## print warning message if profile doesn't meet criteria
            if score < score_cutoff and eValue > eValue_cutoff:
                print '\nWARNING: Bad scores! Profile NOT used.'
                print '\t Name: %s - Score: %f E-value: %f ' \
                      %( hmmName, score, eValue )

            ## add profile to list of profiles to use
            else:
                seqStart = int(hits[h][2])
                seqEnd = int(hits[h][3])
                if hmmHits.has_key(hmmName):
                    hmmHits[hmmName] =  hmmHits[hmmName] +[[seqStart, seqEnd ]]

                else:
                    hmmHits[hmmName] = [ [ seqStart, seqEnd ] ]

        print '\nWill use profile(s):\n ', str( hmmHits ), '\n'

        return hmmHits


    def getHmmProfile( self, hmmName ):
        """
        Get the hmm profile with name hmmName from hmm databse.
        Extract some information about the profile as well as the
        match state emmission scores. Keys of the returned dictionary::
          'AA', 'name', 'NrSeq', 'emmScore', 'accession',
          'maxAllScale', 'seqNr', 'profLength', 'ent', 'absSum'
          
        @param hmmName: hmm profile name
        @type  hmmName: str
        
        @return: dictionary with warious information about the profile
        @rtype: dict
        """
        profileDic = {}
        ## get hmm profile
        try:
            out = os.popen( hmmfetchExe + ' ' + self.hmmdb \
                            + ' ' + hmmName ).read()
        except:
            print 'ERROR getting profile ' + hmmName

        ## write hmm profile to disc (to be used by hmmalign)
        fName = self.hmmFile
        hmm = open( fName, 'w' )
        hmm.write(out)
        hmm.close()

        ## collect some data about the hmm profile
        profileDic['name'] =   hmmName 
        profileDic['profLength'] = \
                  int( string.split(re.findall('LENG\s+[0-9]+', out)[0])[1] )
        profileDic['accession'] = \
                  string.split(re.findall('ACC\s+PF[0-9]+', out)[0])[1] 
        profileDic['NrSeq'] = \
                  int( string.split(re.findall('NSEQ\s+[0-9]+', out)[0])[1] )
        profileDic['AA'] = \
              string.split(re.findall('HMM[ ]+' + '[A-Y][ ]+'*20, out)[0] )[1:]

        ## collect null emmission scores
        pattern = 'NULE[ ]+' + '[-0-9]+[ ]+'*20
        nullEmm = [ float(j) for j in string.split(re.findall(pattern, out)[0])[1:] ]

        ## get emmision scores
        prob=[]
        for i in range(1, profileDic['profLength']+1):
            pattern = "[ ]+%i"%i + "[ ]+[-0-9]+"*20
            e = [ float(j) for j in string.split(re.findall(pattern, out)[0]) ]
            prob += [ e ]

        profileDic['seqNr'] = N.transpose( N.take( prob, (0,),1 ) )
        profileDic['emmScore'] = N.array(prob)[:,1:]

        ## calculate emission probablitities
        emmProb, nullProb = hmmEmm2Prob( nullEmm, profileDic['emmScore'])

        ent = [ N.resize( entropy(e, nullProb), (1,20) )[0] for e in emmProb ]
        profileDic['ent'] = N.array(ent)


        ###### TEST #####

        proba = N.array(prob)[:,1:]

##         # test set all to max score
##         p = proba
##         p1 = []
##         for i in range( len(p) ):
##             p1 += [ N.resize( p[i][N.argmax( N.array( p[i] ) )] , N.shape( p[i] ) ) ]
##         profileDic['maxAll'] = p1

        # test set all to N.sum( abs( probabilities ) )
        p = proba
        p2 = []
        for i in range( len(p) ) :
            p2 += [ N.resize( N.sum( N.absolute( p[i] )), N.shape( p[i] ) ) ]
        profileDic['absSum'] = p2

        # set all to normalized max score 
        p = proba
        p4 = []
        for i in range( len(p) ) :
            p_scale = (p[i] - N.average(p[i]) )/ math.SD(p[i])
            p4 += [ N.resize( p_scale[N.argmax( N.array(p_scale) )] ,
                              N.shape( p[i] ) ) ]
        profileDic['maxAllScale'] = p4

        return profileDic


    def align( self, model, hits ):
        """
        Performs alignment
        If there is more than one hit with the profile, the sequence will
        be subdevided and the alignment will be performed on each part.
        a final merger profile for the profile will be returned.
        
        @param model: model
        @type  model: PDBModel       
        @param hits: list with matching sections from L{searchHmmdb}
        @type  hits: [[int,int]]
                    
        @return: fastaSeq hmmSeq repete hmmGap::
                   fastaSeq - sequence
                   hmmSeq - matching positions in profile
                   repete - number of repetes of the profile
                   hmmGap - list with gaps (deletions in search sequence) for
                            each repete
        @rtype: str, str, int, [int]
        """
        fastaSeq, self.fastaID = self.fasta( model )
        fastaSeq = fastaSeq.split()[1]
        l = len(fastaSeq)

        ## sub sequence alignment
        j = 0
        repete = 0
        hmmGap = []

        for h in hits:
            start = h[0]-1
            stop = h[1]
            sub_fastaSeq, self.fastaID = self.fasta( model, start, stop )

            ## write sub-fasta sequence file
            fName = self.sub_fastaFile
            sub_seq = open( fName, 'w' )
            sub_seq.write( sub_fastaSeq )
            sub_seq.close()

            ## get sub-alignmnet
            sub_fastaSeq, sub_hmmSeq = self.subAlign( self.sub_fastaFile)

            ## remove position scorresponding to insertions in search sequence
            sub_fastaSeq, sub_hmmSeq, del_hmm  = \
                          self.removeGapInSeq( sub_fastaSeq, sub_hmmSeq )
            hmmGap += [ del_hmm ]

            ## rebuild full lenght hmmSeq from sub_hmmSeq
            sub_hmmSeq = '.'*start + sub_hmmSeq + '.'*(l-stop)

            if j == 0:
                hmmSeq = sub_hmmSeq
            if j != 0:
                hmmSeq  = self.mergeHmmSeq( hmmSeq, sub_hmmSeq )

            j+= 1
            repete += 1

        return fastaSeq, hmmSeq, repete, hmmGap


    def removeGapInSeq( self, fasta, hmm ):
        """
        Removes position scorresponding to insertions in search sequence

        @param fasta: search sequence 
        @type  fasta: str
        @param hmm: sequence, matching hmm positions
        @type  hmm: str
                    
        @return: search sequence and profile match sequence with insertions
                 and deletions removed and a list with the deleted positions
        @rtype: str, str, [int]
        """
        new_hmm = []
        new_fasta = []
        del_pos = []

        j = 0
        for i in range( len(fasta) ):
            if string.upper( fasta[i] ) in molUtils.allAA():
                new_hmm += hmm[i]
                new_fasta += fasta[i]
            ## if there is an insertion
            if fasta[i] == '-':
                del_pos += [j]
            ## if there is a deletion
            if hmm[i] == '.':
                j -= 1
            j += 1

        return N.sum(new_fasta), N.sum(new_hmm), del_pos


    def mergeHmmSeq( self, seq1, seq2 ):
        """
        Merges two sequence files into one.
        Multilple hits with one profile cannot overlap!! Overlap == ERROR

        @param seq1: sequence
        @type  seq1: str
        @param seq2: sequence
        @type  seq2: str

        @return: merged sequence or None
        @rtype: str OR None 
        """
        if len(seq1) != len(seq2):
            print 'ERR in mergeHmmSeq:'
            print '\tSequences of different lengths cannot be merged'
            return None
        else:
            result = ''
            for i in range( len(seq1) ):
                ## no match in either
                if seq1[i] == seq2[i] == '.':    
                    result += '.'
                ## match in seq1
                if seq1[i] > seq2[i]:
                    result += seq1[i]
                ## match in seq2
                if seq1[i] < seq2[i]:
                    result += seq2[i]

            return result


    def subAlign( self, file ):
        """
        Align fasta formated sequence to hmm profile.

        @param file: path to hmmalign output file
        @type  file: str
        
        @return: alignment and  matching hmm positions with gaps
        @rtype: str, str
        """
        ## align sequence to hmm profile
        out = os.popen( hmmalignExe + ' -q ' + self.hmmFile + \
                        ' ' + file ).read()

        ## extract search sequence
        fastaSeq = re.findall( self.fastaID + '[ ]+[-a-yA-Y]+', out )
        fastaSeq = string.join([ string.split(i)[1] for i in fastaSeq ], '')

        ## extract hmm sequence
        hmmSeq = re.findall( '#=[A-Z]{2}\s[A-Z]{2}\s+[.x]+', out )
        hmmSeq = string.join([ string.strip( string.split(i)[2] ) for i in hmmSeq ], '')

        return fastaSeq, hmmSeq


    def castHmmDic( self, hmmDic, repete, hmmGap, key ):
        """
        Blow up hmmDic to the number of repetes of the profile used.
        Correct scores for possible deletions in the search sequence.

        @param hmmDic: dictionary from L{getHmmProfile}
        @type  hmmDic: dict
        @param repete: repete information from L{align}
        @type  repete: int
        @param hmmGap: information about gaps from L{align}
        @type  hmmGap: [int]
        @param key: name of scoring method to adjust for gaps and repetes
        @type  key: str
        
        @return: dictionary with information about the profile
        @rtype: dict        
        """
        s = hmmDic[key]

        for i in range( repete ):
            mask = N.ones( len(s) )
            N.put( mask, hmmGap[i], 0 )
            if i == 0:
                score = N.compress( mask, s, 0 )
            if i > 0:
                score = N.concatenate( ( N.compress( mask, s, 0 ), score ) )

        hmmDic[key] = score

        return hmmDic


    def matchScore( self, fastaSeq, hmmSeq, profileDic, key ):
        """
        Get match emmision score for residues in search sequence

        @param fastaSeq: search sequence 
        @type  fastaSeq: str
        @param hmmSeq: sequence, matching hmm positions
        @type  hmmSeq: str
        @param profileDic: from L{castHmmDic}
        @type  profileDic:  dict 
        @param key: name of scoring method
        @type  key: str
        
        @return: list of emmision scores for sequence
        @rtype: [float]
        """
        cons = []
        hmmShift = 0

        for i in range( len( fastaSeq ) ):

            ## if aa is not in hmm profile
            if hmmSeq[i] == '.':
                cons += [0]
                hmmShift += 1

            ## extract emmission value for residue in position i
            if fastaSeq[i] in profileDic['AA'] and  hmmSeq[i] != '.':
                j = profileDic['AA'].index( fastaSeq[i] )
                cons += [ profileDic[key][i - hmmShift][j] ]

        return cons


    def __score( self, model, key, hmmNames=None ):
        """
        Get match emmission scores for search sequence.
        If profile name(s) are provided, no search will be performed.

         - If names and positions of hmm profiles is B{NOT} provided
           B{search performed} -> score (array), hmmNames (dictionary)
           
         - If names and positions of hmm profiles is provided
           B{NO search performed} -> score (array), hmmNames
           (dictionary - same as input)

        @param model: model
        @type  model: PDBModel
        @param key: name of scoring method
        @type  key: str    
        @param hmmNames: profile name OR None (default: None)
        @type  hmmNames: str OR None
        
        @return: score, hmmNames
        @rtype: array, dict
        """
        if not hmmNames:
            ## get profile for model
            searchResult, searchHits = self.searchHmmdb( model )
            hmmNames = self.selectMatches( searchResult, searchHits )
        else:
            searchResult = self.searchHmmdb( model, noSearch=1 )
            hmmNames = hmmNames

        ## retrieve hmm model(s)
        result = None

        for name in hmmNames.keys():
            ## retrieve hmm model
            hmmDic = self.getHmmProfile( name )
            print 'Hmm profile ' + str(name) + ' with accession number ' \
                  + str(hmmDic['accession']) + ' retrieved'

            ## align sequence with model
            fastaSeq, hmmSeq, repete, hmmGap = self.align( model,
                                                           hmmNames[ name ] )
            ## cast hmm model
            hmmDic = self.castHmmDic( hmmDic, repete, hmmGap, key )

            ## Hmmer profile match scores for sequence
            cons = self.matchScore( fastaSeq, hmmSeq, hmmDic, key )

            if result:
                result = self.mergeProfiles( result, cons )
            else:
                result = cons

        return result, hmmNames


    def scoreAbsSum( self, model, hmmNames=None ):
        return self.__score( model, key='absSum', hmmNames=hmmNames )


    def scoreMaxAll( self, model, hmmNames=None ):
        return self.__score( model, key='maxAllScale', hmmNames=hmmNames )


    def scoreEntropy( self, model, hmmNames=None ):
        return self.__score( model, key='ent', hmmNames=hmmNames )


    def __list2array( self, lstOrAr ):
        if type( lstOrAr ) == list:
            return N.array( lstOrAr )
        return lstOrAr


    def mergeProfiles( self, p0, p1, maxOverlap=3 ):
        """
        Merge profile p0 with profile p1, as long as they overlap in
        at most maxOverlap positions

        @param p0: profile
        @type  p0: [float]
        @param p1: profile
        @type  p1: [float]
        @param maxOverlap: maximal allowed overlap between profiles
        @type  maxOverlap: int
        
        @return: array
        @rtype: 
        """
        p0 = self.__list2array( p0 )
        p1 = self.__list2array( p1 )

        overlap = N.greater( N.greater(p0,0) + N.greater(p1,0), 1 )

        if N.sum( overlap ) <= maxOverlap:
            ## one of the two profiles will in most cases not belong to these
            ## positions. We can't decide which one is wrong, let's eliminate
            ## both values. Alternatively we could keep one, or the average, ..
            N.put( p1, N.nonzero( overlap ), 0 )
            N.put( p0, N.nonzero( overlap ), 0 )

            p0 = p0 + p1

        return p0


    def cleanup( self ):
        """
        remove temp files
        """
        del_lst = [ self.sub_fastaFile, self.hmmFile, self.fastaFile ]
        for d in del_lst:
            try:
                os.remove( d )
            except:
                pass


####################################
## Testing

if __name__ == '__main__':

    import tools as T
    import glob
    from PDBModel import PDBModel

    print "Loading PDB..."
    f = glob.glob( T.testRoot()+'/lig_pcr_00/pcr_00/*_1_*pdb' )[1]
    m = PDBModel(f)
    model = m.compress( m.maskProtein() )

    ## initiate and check database status
    a = Hmmer( hmmdb = settings.hmm_db )
    a.checkHmmdbIndex()

    ## scoring methods to use
    method = [ 'emmScore', 'ent', 'maxAll', 'absSum', 'maxAllScale' ]

    ## search
    searchMatches, searchHits = a.searchHmmdb( model )
    hmmNames = a.selectMatches( searchMatches, searchHits )

    cons = []
    result = None

    for name in hmmNames.keys():

        ## retrieve hmm model
        hmmDic = a.getHmmProfile( name )

        ## align sequence with model
        fastaSeq, hmmSeq, repete, hmmGap = a.align( model, hmmNames[ name ] )

        ## cast hmm model
        hmmDic_cast = a.castHmmDic( hmmDic, repete, hmmGap, method[0] )

        ## Hmmer profile match scores for sequence
        cons = a.matchScore( fastaSeq, hmmSeq, hmmDic_cast, method[0] )

        ## If there are more than one profile in the model, merge to one. 
        if result:
            result = a.mergeProfiles( result, cons )
        else:
            result = cons

    ## cleanup
    a.cleanup()
