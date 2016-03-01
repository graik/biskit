## Automatically adapted for numpy.oldnumeric Mar 26, 2007 by alter_code1.py

## Class conservation
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
## $Revision$
## last $Author$
## last $Date$

"""
Search Hmmer Pfam database and retrieve conservation data.
"""

import tempfile
import re, string
import types, os.path
import numpy as N
import molUtils
import settings

import Biskit.mathUtils as math
from Biskit.Errors import BiskitError
import Biskit.tools as T
import Biskit.molTools as MT
import Biskit.molUtils as MU
from Biskit import Executor, TemplateError, PDBModel, StdLog, EHandler


## executables
## hmmpfamExe = settings.hmmpfam_bin
## hmmfetchExe = settings.hmmfetch_bin
hmmalignExe = settings.hmmalign_bin
## hmmindexExe = settings.hmmindex_bin
hmmDatabase = settings.hmm_db


class HmmerError( BiskitError ):
    pass


class HmmerdbIndex( Executor ):
    """
    Checks if the hmm database has been indexed or not,
    if not indexing will be done.
    """
    
    def __init__( self, hmmdb, **kw ):
        """
        @param hmmdb: Pfam hmm database
        @type  hmmdb: str
        """
        Executor.__init__( self, 'hmmindex', args='%s'%hmmdb, **kw )

        if not os.path.exists(hmmdb+'.ssi'):
            if self.verbose:
                self.log.writeln(
                    'HMMINDEX: Indexing hmm database. This will take a while')
            
            self.run()


class HmmerSearch( Executor ):
    """
    Search hmm database (using the hmmpfam program) with a sequence
    in fasta format.
    If the profile names have been provided - skip the search and
    only write the temporary sequence files.
    """

    def __init__( self, target, hmmdb=settings.hmm_db, noSearch=None, **kw ):
        """
        @param target: fasta sequence, fasta file, or PDBModel 
        @type  target: PDBModel or str (fasta file) or [ str ] (fasta lines)
        @param hmmdb: Pfam hmm database
        @type  hmmdb: str
        @param noSearch: don't perform a seach
        @type  noSearch: 1 OR None
        """
        self.hmmdb = hmmdb

        Executor.__init__( self, 'hmmpfam', f_in= tempfile.mktemp('.fasta'),
                           catch_out=1, **kw )

        self.target = target
        
        self.fastaID = ''
        
        if noSearch:
            if self.verbose:
                self.log.writeln(
                    'Profiles provided - No search will be performed.')


    def __verify_fasta(self, target ):
        """
        Verify that a given file or string is in Fasta format.
        The definition used for a fasta file here is that:
         - first line starts with '>'
         - the following sequence lines are not longer that 80 characters
         - the characters has to belong to the standard amino acid codes

        @param target: name of fasta file OR file contents as list of strings
        @type  target: str OR [str]

        @return: conforms to the fsata format
        @rtype: True/False
        """
        if type(target) != types.ListType:
            if os.path.exists( target ):
                f = open( target, 'r' )
                target = f.readlines()

        if target[0][0] != '>':
            self.log.writeln('Fasta format does not contain description line.')
            return False

        for i in range( 1, len(target) ):
            if len( target[i] ) >= 80:
                self.log.writeln(
                    'Fasta sequence lines longer that 80 characters')
                return False

            for j in target[i]:
                aa_codes = MU.aaDicStandard.values() + [ '\n' ]
                if not j.upper() in aa_codes:
                    self.log.writeln('Invalid amino acid code: %s'%j.upper())
                    return False

        return True

           
    def prepare( self ):
        """
        If the target is a PDBModel its sequence will be written
        to disc as a fasta file.
        If it is already fasta file the path to the file will be passed on.
        """
        try:
            ## if target is a PDBModel or simple sequence string
            if isinstance(self.target, PDBModel):
                fastaSeq, self.fastaID = MT.fasta( self.target )
                ## write fasta sequence file
                seq = open( self.f_in, 'w' )
                seq.write( fastaSeq )
                seq.close()
            
            elif self.__verify_fasta(self.target):

                if type(self.target) is list:
                    seq = open( self.f_in, 'w' )
                    seq.writelines( self.target )
                    seq.close()
                else:
                    ## else assume it is a fasta sequence file   
                    self.f_in = self.target
            
            else:
                raise HmmerError, 'cannot interpret target %r' % self.target

        except OSError, why:
            msg = 'Cannot write temporary fasta file %s' % self.f_in
            msg += '\nError: %r' % why
            raise HmmerError, msg
        
        ## fix bug 2816430
        self.args = ' %s %s' % (self.hmmdb, self.f_in)
        

    def parse_result( self ):
        """
        Parse the output from hmmpfam.
        
        @return: dictionary with profile names as keys and a list of
                 lists containing information about the range where the
                 profile matches the sequence
        @rtype: dict, [list]
        """
        matches = {}
        hits = []

        ## check that the outfut file is there and seems valid
        if not os.path.exists( self.f_out ):
            raise HmmerError,\
                  'Hmmersearch result file %s does not exist.'%self.f_out
        
        if T.fileLength( self.f_out ) < 10:
            raise HmmerError,\
                  'Hmmersearch result file %s seems incomplete.'%self.f_out 

        try:
            lines = open( self.f_out, 'r' ).readlines()
            while lines:
                l = lines.pop(0)
                ## get names and descriptions of matching profiles
                if re.match('^-{8}\s{7,8}-{11}.+-{3}$', l):
                    m = string.split( lines.pop(0) )
                    while len(m) != 0:
                        matches[m[0]] =  m[1:] 
                        m = string.split( lines.pop(0) )

                ## get hits, scores and alignment positions
                if re.match('^-{8}\s{7,8}-{7}\s-{5}\s-{5}.+-{7}$', l):
                    h = string.split( lines.pop(0) )
                    while len(h) != 0:
                        hits += [ h ] 
                        h = string.split( lines.pop(0) )
                    break
                
        except Exception, why:
            raise HmmerError,\
                  'ERROR parsing hmmpfam search result: %s'%self.f_out +\
                  '\n%r' % why
        
        return matches, hits
        

    def finish( self ):
        """
        Overrides Executor method
        """
        Executor.finish( self )
        self.result = self.parse_result( )



class HmmerProfile( Executor ):
    """
    Get the hmm profile with name hmmName from hmm databse
    (using the program hmmfetch) and calculate 3 different conservation
    scores.

    Use::

      p = HmmerProfile( 'FH2' )
      profileDic = p.run()

    Will fetch the Pfam profile 'FH2'.

    profileDic then has the following entries:

      * 'name'
      * 'accession'
      * 'NrSeq'
      * 'profLength' -- profile length
      * 'seqNr' -- ?? list 1 .. NrSeq 
      * 'AA' -- amino acids in same order as columns in 'emmScore'
      * 'emmScore' -- (profLength x 20) array of emmission scores
      * 'ent' -- Kulbach-Leibler distance conservation measure (best measure)
      * 'absSum' -- alternative conservation measure
      * 'maxAllScale' -- alternative conservation measure

    Use 
    >>> p = HmmerProfile( 'FH2', f_out='your_new.hmm' )
    to capture the HMM file, which is otherwise cleaned up.
    """
    
    def __init__( self, hmmName, hmmdb=settings.hmm_db, **kw ):
        """
        @param hmmName: hmm profile name
        @type  hmmName: str
        @param hmmdb: Pfam hmm database
        @type  hmmdb: str
        **kw - all Executor parameters, in particular:
        @param f_out: target file name for profile.hmm, prevents its deletion
        @type  f_out: str
        @param debug:
        @param log:
        @param ...
        """
        self.hmmName = hmmName

        Executor.__init__( self, 'hmmfetch',
                           args=' %s %s'%(hmmdb, hmmName), **kw )

        
    def parse_result( self ):
        """
        Extract some information about the profile as well as the
        match state emmission scores. Keys of the returned dictionary::
          'AA', 'name', 'NrSeq', 'emmScore', 'accession',
          'maxAllScale', 'seqNr', 'profLength', 'ent', 'absSum'
          
        @return: dictionary with warious information about the profile
        @rtype: dict
        """
        ## check that the outfut file is there and seems valid
        if not os.path.exists( self.f_out ):
            raise HmmerError,\
                  'Hmmerfetch result file %s does not exist.'%self.f_out
        
        if T.fileLength( self.f_out ) < 10:
            raise HmmerError,\
                  'Hmmerfetch result file %s seems incomplete.'%self.f_out
        
        profileDic = {}

        ## read result
        hmm = open( self.f_out, 'r')
        out = hmm.read()
        hmm.close()

        ## collect some data about the hmm profile
        profileDic['name'] =  self.hmmName 
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
        emmProb, nullProb = self.hmmEmm2Prob( nullEmm, profileDic['emmScore'])

        ent = [ N.resize( self.entropy(e, nullProb), (1,20) )[0] for e in emmProb ]
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


    def hmmEmm2Prob( self, nullEmm, emmScore ):
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
        ## see http://www.ebc.ee/WWW/hmmer2-html/node26.html
        emmProb = nullProb * N.power( 2, ( emmScore/1000.0) )

        return emmProb, nullProb


    def entropy( self, emmProb, nullProb ):
        """ 
        Calculate the Kullback-Leibler distance between the observed and the
        background amino acid distribution at a given position. High values mean
        high conservation. Empty (all 0) emmission probabilities yield score 0.
        See also:BMC Bioinformatics. 2006; 7: 385
        
        emmProb & nullProb is shape 1,len(alphabet)

        @param emmProb: emmission probabilities
        @type  emmProb: array
        @param nullProb: null probabilities
        @type  nullProb: array

        @return: relative entropy score
        @rtype:  float
        """
        ## avoid log error
        if N.sum( emmProb ) == 0.:
            return 0.

        return N.sum( emmProb * N.log(emmProb/nullProb) )

            
    def fail( self ):
        """
        Called if external program failed, override! 
        """
        raise HmmerError,\
              'Hmmerfetch failed retrieving profile: %s'%self.hmmName

        
    def finish( self ):
        """
        Overrides Executor method
        """
        Executor.finish( self )
        self.result = self.parse_result( )


class HmmerAlign( Executor ):
    """
    Align fasta formated sequence to hmm profile (using hmmalign).
    """
    
    def __init__( self, hmmFile, fastaFile, fastaID, **kw ):
        """
        @param hmmFile: path to hmm file (profile)
        @type  hmmFile: str
        @param fastaFile: path to fasta search sequence
        @type  fastaFile: str
        @param fastaID: fasta id of search sequence
        @type  fastaID: str     
        """
        self.fastaID = fastaID
        self.hmmFile = hmmFile
        self.fastaFile = fastaFile
        
        assert T.fileLength( self.hmmFile ) > 10, \
               'input HMM file missing or empty'
        
        Executor.__init__( self, 'hmmalign',
                           args=' -q %s %s'%(hmmFile, fastaFile), **kw )

        
    def parse_result( self ):
        """
        Align fasta formated sequence to hmm profile.
        
        @return: alignment and  matching hmm positions with gaps
        @rtype: str, str
        """
        ## check that the outfut file is there and seems valid
        if not os.path.exists( self.f_out ):
            raise HmmerError,\
                  'Hmmeralign result file %s does not exist.'%self.f_out
        
        if T.fileLength( self.f_out ) < 1:
            raise HmmerError,\
                  'Hmmeralign result file %s seems incomplete.'%self.f_out
        
        ## read result
        hmm = open( self.f_out, 'r')
        out = hmm.read()
        hmm.close()
        
        ## extract search sequence
        fastaSeq = re.findall( self.fastaID + '[ ]+[-a-yA-Y]+', out )
        fastaSeq = string.join([ string.split(i)[1] for i in fastaSeq ], '')

        ## extract hmm sequence
        hmmSeq = re.findall( '#=[A-Z]{2}\s[A-Z]{2}\s+[.x]+', out )
        hmmSeq = string.join([ string.strip( string.split(i)[2] ) for i in hmmSeq ], '')

        return fastaSeq, hmmSeq


    def fail( self ):
        """
        Called if external program failed, override! 
        """
        raise HmmerError,\
              'hmmeralign failed aligning sequence file %s with profile %s'\
              %(self.fastaFile ,self.hmmFile)
    

    def finish( self ):
        """
        Overrides Executor method
        """
        Executor.finish( self )
        self.result = self.parse_result( )

        
class Hmmer:
    """
    Search Hmmer Pfam database and retrieve conservation score for model
    """

    def __init__(self, hmmdb=hmmDatabase, verbose=1, log=StdLog(),
                 debug=False ): 
        """
        @param hmmdb: Pfam hmm database
        @type  hmmdb: str
        @param verbose: verbosity level (default: 1)
        @type  verbose: 1|0
        @param log: Log file for messages [STDOUT]
        @type  log: Biskit.LogFile
        @param debug: don't cleanup temporary files [False]
        @type  debug: 1|0
        
        """
        self.hmmdb = hmmdb
        self.verbose = verbose
        self.log = log or StdLog()
        self.tempDir = settings.tempDirShared
        self.debug = debug

        self.fastaID = ''

        self.hmmFile = tempfile.mktemp('.hmm', dir=self.tempDir)
        self.fastaFile = tempfile.mktemp('.fasta', dir=self.tempDir)
        self.sub_fastaFile = tempfile.mktemp('_sub.fasta', dir=self.tempDir)


    def checkHmmdbIndex( self ):
        """
        Checks if the hmm database has been indexed or not,
        if not indexing will be done.
        """
        idx = HmmerdbIndex( self.hmmdb, verbose=self.verbose, log=self.log )


    def searchHmmdb( self, target, noSearch=None ):
        """
        Search hmm database with a sequence in fasta format.
        If the profile names have been provided - skip the search and
        only write the temporary sequence files.
          
        @param target: sequence file or PDBModel
        @type  target: PDBModel or fasta file
        @param noSearch: don't perform a seach
        @type  noSearch: 1 OR None
         
        @return: dictionary witn profile names as keys and a list of
                 lists containing information about the range where the
                 profile matches the sequence
        @rtype: dict, [list]
        """
        if self.verbose:
            self.log.writeln(
                '\nSearching hmm database, this will take a while...')

        search = HmmerSearch( target, self.hmmdb, verbose=self.verbose,
                              log=self.log, debug=self.debug )
        matches, hits = search.run()
        
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
            if self.verbose: self.log.writeln('\nNo matching profiles found')
            return None
        if len(matches) == 1:
            if self.verbose:
                self.log.writeln('\nFound only one matching profile.')
        if len(matches) > 1:
            if self.verbose: self.log.writeln('\nFound more than one hit.\n')

        try:
            if self.verbose:
                self.log.writeln(
                    '\tName -- Nr hits -- Score -- E-value -- Description')
                for k in matches.keys():
                    self.log.writeln( '\t'+k+' --'+matches[k][-1] +' -- '+\
                          matches[k][-3]+ ' -- '+ matches[k][-2]+' -- '+\
                          string.join(matches[k][0:-3]))
        except:
            pass

        ## select to use profiles from searchHits  given the given criteria
        hmmHits = {}

        for h in range( len(hits) ):

            hmmName = hits[h][0]
            score = float( hits[h][-2] )
            eValue = float( hits[h][-1] )

            ## print warning message if profile doesn't meet criteria
            if score < score_cutoff and eValue > eValue_cutoff \
                   and self.verbose:
                self.log.writeln('\nWARNING: Bad scores! Profile NOT used.')
                self.log.writeln('\t Name: %s - Score: %f E-value: %f ' \
                      %( hmmName, score, eValue ) )

            ## add profile to list of profiles to use
            else:
                seqStart = int(hits[h][2])
                seqEnd = int(hits[h][3])
                if hmmHits.has_key(hmmName):
                    hmmHits[hmmName] =  hmmHits[hmmName] +[[seqStart, seqEnd ]]

                else:
                    hmmHits[hmmName] = [ [ seqStart, seqEnd ] ]

        if self.verbose:
            self.log.writeln( '\nWill use profile(s):\n '+str( hmmHits ))

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
        profile = HmmerProfile( hmmName, self.hmmdb,
                                # giving explicit f_out prevents deletion
                                f_out=self.hmmFile,
                                verbose=self.verbose,
                                log=self.log,
                                debug=self.debug )

        r = profile.run()
        assert os.path.exists( self.hmmFile )
        return r

    def align( self, model, hits ):
        """
        Performs alignment
        If there is more than one hit with the profile, the sequence will
        be subdevided and the alignment will be performed on each part.
        a final merger profile for the profile will be returned.
        
        @param model: structure model
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
        try:

            fastaSeq, self.fastaID = MT.fasta( model )
            fastaSeq = fastaSeq.split()[1]
            l = len(fastaSeq)

            ## sub sequence alignment
            j = 0
            repete = 0
            hmmGap = []

            for h in hits:
                start = h[0]-1
                stop = h[1]
                sub_fastaSeq, self.fastaID = MT.fasta( model, start, stop )

                ## write sub-fasta sequence file
                fName = self.sub_fastaFile
                sub_seq = open( fName, 'w' )
                sub_seq.write( sub_fastaSeq )
                sub_seq.close()

                ## get sub-alignmnet
                align = HmmerAlign( self.hmmFile, self.sub_fastaFile,
                                    self.fastaID, verbose=self.verbose,
                                    log=self.log)

                sub_fastaSeq, sub_hmmSeq = align.run()

                ## remove positions corresponding to insertions in
                ## search sequence
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

        except IOError, e:
            raise HmmerError, "Error creating temporary file %r. "%e.filename+\
                  "Check that Biskit.settings.temDirShared is accessible!\n"+\
                  "See also .biskit/settings.cfg!"


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

        return ''.join(new_fasta), ''.join(new_hmm), del_pos


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
            EHandler.warning( 'ERR in mergeHmmSeq:\n' +\
                         '\tSequences of different lengths cannot be merged')
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
##         out = os.popen( hmmalignExe + ' -q ' + self.hmmFile + \
##                         ' ' + file ).read()
        out = os.popen( self.exe.bin + ' -q ' + self.hmmFile + \
                        ' ' + file ).read()

        ## extract search sequence
        fastaSeq = re.findall( self.fastaID + '[ ]+[-a-yA-Y]+', out )
        fastaSeq = string.join([ string.split(i)[1] for i in fastaSeq ], '')

        ## extract hmm sequence
        hmmSeq = re.findall( '#=[A-Z]{2}\s[A-Z]{2}\s+[.x]+', out )
        hmmSeq = string.join(
            [ string.strip( string.split(i)[2] ) for i in hmmSeq ], '')

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

        assert len(hmmSeq) == len( fastaSeq ), \
               'Length of HMM profile does not match length of sequence.\n' + \
               'Check for unusual residues in input model!'

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
        #else:
            #searchResult = self.searchHmmdb( model, noSearch=1 )
            #hmmNames = hmmNames

        ## retrieve hmm model(s)
        result = None

        for name in hmmNames.keys():
            ## retrieve hmm model
            hmmDic = self.getHmmProfile( name )
            if self.verbose:
                self.log.writeln('Hmm profile ' + str(name) +\
                                 ' with accession number ' \
                                 + str(hmmDic['accession']) + ' retrieved')

            ## align sequence with model
            fastaSeq, hmmSeq, repete, hmmGap = self.align( model,
                                                           hmmNames[ name ] )
            ## cast hmm model
            hmmDic = self.castHmmDic( hmmDic, repete, hmmGap, key )

            ## Hmmer profile match scores for sequence
            cons = self.matchScore( fastaSeq, hmmSeq, hmmDic, key )

            if result is not None:
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
            N.put( p1, N.nonzero( overlap )[0], 0 )
            N.put( p0, N.nonzero( overlap )[0], 0 )

            p0 = p0 + p1

        return p0


    def cleanup( self ):
        """
        remove temp files
        """
        if not self.debug:
            del_lst = [ self.sub_fastaFile, self.hmmFile, self.fastaFile ]
            for d in del_lst:
                T.tryRemove( d )


#############
##  TESTING        
#############
import Biskit.test as BT
        
class Test(BT.BiskitTest):
    """Hmmer test"""

    TAGS = [ BT.EXE ]
    
    M = None

    def prepare( self ):
        from Biskit import PDBModel
        import Biskit.tools as T

        if not self.M:
            if self.local: print "Loading PDB...",
            self.M = PDBModel( T.testRoot()+'/lig/1A19.pdb')
            self.M = self.M.compress( self.M.maskProtein() )
            if self.local: print "Done"
    
    def test_HmmSearch( self ):
        """HmmSearch test"""
        self.searcher = HmmerSearch( self.M, 
                                   verbose=self.local,
                                   log=self.log )
        self.matches, self.hits = self.searcher.run()

    def test_HmmProfile( self ):
        """HmmerProfile test"""
        profile = HmmerProfile( 'FH2', verbose=self.local, log=self.log)
        self.profileDic = profile.run()

    def test_HmmerFasta( self ):
        """Hmmer test (search from fasta)"""
        h = Hmmer(hmmdb=settings.hmm_db)
        h.checkHmmdbIndex()

        self.searchMatches, self.searchHits = h.searchHmmdb(
            T.testRoot()+'/Mod/project/target.fasta' )
        
        self.assertTrue( len( self.searchHits ) > 3 )

    def test_HmmerModel( self):
        """Hmmer test (search from model)"""
        
        ## initiate and check database status
        self.hmmer = Hmmer( hmmdb=settings.hmm_db, verbose=self.local,
                            log=self.log )
        self.hmmer.checkHmmdbIndex()

        ## scoring methods to use
        method = [ 'emmScore', 'ent', 'maxAll', 'absSum', 'maxAllScale' ]

        ## search
        searchMatches, searchHits = self.hmmer.searchHmmdb( self.M )
        hmmNames = self.hmmer.selectMatches( searchMatches, searchHits )
  ##      hmmNames = {'Barstar': [[1, 89]]}
        
        self.cons = []
        self.result = None

        for name in hmmNames.keys():

            ## retrieve hmm model
            hmmDic = self.hmmer.getHmmProfile( name )

            ## align sequence with model
            fastaSeq, hmmSeq, repete, hmmGap = \
                  self.hmmer.align( self.M, hmmNames[ name ] )

            ## cast hmm model
            hmmDic_cast = \
                  self.hmmer.castHmmDic( hmmDic, repete, hmmGap, method[0] )

            ## Hmmer profile match scores for sequence
            self.cons = self.hmmer.matchScore( fastaSeq, hmmSeq,
                                          hmmDic_cast, method[0] )

            ## If there are more than one profile in the model, merge to one.
            ## Note: this can be problematic as scores depend on family size
            if self.result:
                self.result = self.hmmer.mergeProfiles( self.result, self.cons )
            else:
                self.result = self.cons

        self.hmmer.cleanup()

        self.assertEqual( self.result, self.EXPECTED )


    #: Hmmer emission scores
    EXPECTED = [2581.0, 3583.0, 1804.0, 2596.0, 3474.0, 2699.0, 3650.0, 2087.0, 2729.0, 2450.0, 2412.0, 2041.0, 3474.0, 1861.0, 2342.0, 2976.0, 5124.0, 2729.0, 2202.0, 2976.0, 3583.0, 2202.0, 2103.0, 2976.0, 1922.0, 2132.0, 4122.0, 2403.0, 4561.0, 4561.0, 3650.0, 2087.0, 4001.0, 2976.0, 3860.0, 3260.0, 2976.0, 6081.0, 3860.0, 5611.0, 2976.0, 3609.0, 3650.0, 6081.0, 3343.0, 2403.0, 3288.0, 4122.0, 2976.0, 2322.0, 2976.0, 1995.0, 4378.0, 2706.0, 2665.0, 4186.0, 3539.0, 2692.0, 3270.0, 2302.0, 2604.0, 2132.0, 2118.0, 2380.0, 2614.0, 2170.0, 3260.0, 2403.0, 1964.0, 3343.0, 2976.0, 2643.0, 3343.0, 2714.0, 2591.0, 3539.0, 3260.0, 2410.0, 1809.0, 3539.0, 2111.0, -774.0, 3860.0, 2450.0, 2063.0, 3474.0, 3474.0, 2057.0, 1861.0]


if __name__ == '__main__':

    BT.localTest()
    
    



