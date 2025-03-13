## numpy-oldnumeric calls replaced by custom script; 09/06/2016
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
collect and manage information about a docking result
"""

from copy import deepcopy, copy
from numpy import ndarray as arraytype
from difflib import SequenceMatcher

import biskit
import biskit.core.oldnumeric as N0
from biskit import PDBModel, molUtils, mathUtils, StdLog, EHandler
import biskit.tools as t

from biskit import PDBDope
from biskit.exe import Prosa2003

class ComplexError(Exception):
    pass


class Complex:
    """
    Manage receptor and ligand structure, store additional information,
    calculate some properties.
    """

    def __init__(self, rec_model=None,lig_model=None, ligMatrix=None,info={} ):
        """
        @param rec_model: model of original receptor conformation
        @type  rec_model: PDBModel OR XplorModel
        @param lig_model: model of original ligand conformation
        @type  lig_model: PDBModel OR XplorModel
        @param ligMatrix: Numeric array 4 by 4, ligand transformation matrix
        @type  ligMatrix: matrix
        @param info: optional dictionary with additional infos
                     e.g. {'eshape':-123.3, 'rms':12.2 }
        @type  info: dict
        """
        self.rec_model = rec_model  # XplorModel object for receptor
        self.lig_model = lig_model  #    "            for ligand
        self.lig_transformed = None #    "     with transformed coordinates
        self.pw_dist = None         # cached pw atom distances rec x lig

        self.info = { 'date':t.dateSortString() } ## default info record
        self.info.update( info )

        self.ligandMatrix = ligMatrix
        if self.ligandMatrix is None:
            self.ligandMatrix = N0.array([ [1,  0,  0, 0],
                                        [0,  1,  0, 0],
                                        [0,  0,  1, 0],
                                        [0,  0,  0, 1],], N0.Float32)

        ## compressed by slim
        self.contacts = None

        ## version as of creation of this object
        self.initVersion = biskit.__version__


    def __getstate__( self ):
        """
        Called before pickling.
        """
        self.slim()
        return self.__dict__


    def __getitem__( self, key ):
        """
        @return: C.info[ key ]
        @rtype: dict
        """
        return self.info[ key ]


    def __setitem__( self, key, v ):
        self.info[ key ] = v


    def __delitem__( self, key ):
        del self.info[ key ]


    def __contains__( self, key ):
        return key in self.info


    def __setstate__(self, state ):
        """
        called for unpickling the object.
        """
        self.__dict__ = state
        self.ligandMatrix = N0.array( self.ligandMatrix,N0.Float32 )
        ## backwards compability
        self.__defaults() 


    def __defaults(self ):
        """
        backwards compatibility to earlier pickled models
        """
        self.pw_dist = getattr( self, 'pw_dist', None )


    def keys( self ):
        """
        Keys of info dictionary.

        @return: list of keys
        @rtype: [str]
        """
        return list(self.info.keys())


    def has_key( self, k ):
        """
        Check if info dictionary has key.
        
        @param k: key to test
        @type  k: str
        
        @return: dict has key
        @rtype: 1|0
        """
        return k in self.info


    def values( self, keys=[], default=None ):
        """
        Use::
           values( keys=None, default=None ) -> list of info dict values
        
        @param keys: give only values of these keys
        @type  keys: [str]
        @param default: only used with keys!=[], default value for missing keys
        @type  default: any
        
        @return: all values (default) or values of requested keys (ordered)
        @rtype: [any]
        """
        if not keys:
            return list(self.info.values())
        return [ self.get( k, default ) for k in keys ]


    def get( self, key, default=None):
        """
        Use::
           C.get(k[,default]) -> C.info[k] if C.info.has_key(k), else
                                 default defaults to None.

        @param key: key to call
        @type  key: str
        @param default: default value if dictionary doesn't have key
                        (default: None)
        @type  default: any

        @return: dictionary value of key OR else default
        @rtype: any OR None    
        """
        return self.info.get( key, default )


    def rec(self):
        """
        Return XplorModel object.

        @return: receptor model
        @rtype: XplorModel
        """
        return self.rec_model


    def lig(self, force=0, cache=1 ):
        """
        Return ligand structure transformed. Use cached object, if
        available.

        @param force: force new transformation (default: 0)
        @type  force: 1|0
        @param cache: cache transformed model (default: 1)
        @type  cache: 1|0

        @return: ligand model
        @rtype: XplorModel
        """
        try:
            lig = self.lig_transformed
        except:
            lig = None

        if lig is None or force:
            ## get transformation matrix and apply it
            lig = self.lig_model.transform( *self.ligMatrix() )

            if cache:
                self.lig_transformed = lig

        return lig


    def model(self):
        """
        Model with both receptor and ligand.
        
        @return: single PDBModel with first rec then lig
        @rtype: XplorModel
        """
        return self.rec().concat( self.lig() )


    def ligMatrix(self):
        """
        Return transformation matrix to transform original Ligand
        XplorModel into docked conformation.
        
        @return: tuple with rotation and transformation matrix
        @rtype: array, vector
        """
        a = self.ligandMatrix
        ## return tuple of rotation matrix and translation vector
        return (a[0:3,0:3], a[0:3, 3])


    def setLigMatrix( self, m ):
        """
        Set a ligand translation/rotation matrix.
        
        @param m: translation array 4x4 with rotation
        @type  m: array 
        """
        self.ligandMatrix = m
        self.lig_transformed = None


    def getInfo(self):
        """
        Return contents of the info dictionary.

        @return: info dictionary
        @rtype: dict
        """
        return self.info


    def writeLigand(self, fname, ter=0 ):
        """
        Write transformed ligand to pdb file. Remove TER record for
        xplor usage (unless removeTer!=0).
        
        @param fname: output filename
        @type  fname: str
        @param ter: option for how to treat the TER statement::
                      - 0, don't write any TER statements (default)
                      - 1, restore original TER statements
                      - 2, put TER between all detected chains
        @type  ter: 0|1|2
        """
        pdb = self.lig()
        pdb.writePdb( t.absfile( fname ), ter )


    def writeReceptor(self, fname, ter=0 ):
        """
        Write receptor to pdb without TER statement (unless removeTer!=0).
        
        @param fname: output file name
        @type  fname: str
        @param ter: option for how to treat the TER statement::
                      - 0, don't write any TER statements (default)
                      - 1, restore original TER statements
                      - 2, put TER between all detected chains
        @type  ter: 0|1|2
        """
        self.rec().writePdb( t.absfile( fname), ter )


    def writeComplex(self, fname, ter=0):
        """
        Write single pdb containing both receptor and ligand
        separated by END but not TER (unless ter!=0).
        
        @param fname: filename of pdb
        @type  fname: str
        @param ter: option for how to treat the TER statement::
                      - 0, don't write any TER statements (default)
                      - 1, restore original TER statements
                      - 2, put TER between all detected chains
        """
        self.model().writePdb( t.absfile( fname ), ter )


    def slim(self):
        """
        Remove coordinates and atoms of ligand and receptor from memory,
        if they can be restored from file, compress contact matrix.
        @note: CALLED BEFORE PICKLING
        """
        self.lig_transformed = None
        self.pw_dist = None

##         self.ligandMatrix = self.ligandMatrix.tolist()

        if 'matrix' in self.info:
            del self.info['matrix']

        ## compress contact matrix array
        if self.contacts is not None and \
               len(N0.shape( self.contacts['result'] ) )==2:
            m = self.contacts['result']
            self.contacts['shape'] = N0.shape( m )

            self.contacts['result'] = N0.nonzero( N0.ravel( m ) ).astype(N0.Int32)


    def contactsOverlap(self, ref, cutoff=None):
        """
        Fraction of overlapping B{residue-residue} contacts between this and
        reference complex.
        
        @param ref: reference complex
        @type  ref: Complex
        @param cutoff: maximal atom-atom distance, None .. previous setting
        @type  cutoff: float
        
        @return: fraction of contacts shared between this and ref
                 (normalized to number of all contacts)
        @rtype: float
        """
        equal = N0.logical_and(self.resContacts( cutoff=cutoff ),
                            ref.resContacts( cutoff=cutoff ) )
        total = N0.logical_or( self.resContacts(cutoff),
                              ref.resContacts(cutoff) )

        return N0.sum(N0.sum( equal )) * 1.0 / N0.sum(N0.sum( total ))


    def contactsShared(self, reference, cutoff=None):
        """
        Number of equal B{residue-residue} contacts in this and
        reference complex.
        
        @param reference: reference complex
        @type  reference: Complex
        @param cutoff: cutoff for atom-atom contact to be counted
        @type  cutoff: float
        @return: the number or residue-residue contacts that are common to
                 both this and reference::
                   abs( N0.sum( N0.sum( contactMatrix_a - contactMatrix_b )))
        @rtype: int
        """
        equality = N0.logical_and(self.resContacts( cutoff=cutoff ),
                               reference.resContacts( cutoff=cutoff ) )
        return abs(N0.sum(N0.sum( equality )))


    def contactsDiff(self, ref, cutoff=None):
        """
        Number of different B{residue-residue} contacts in this and
        reference complex.
        
        @param ref: to compare this one with
        @type  ref: Complex
        @param cutoff: maximal atom-atom distance, None .. previous setting
        @type  cutoff: float
        
        @return: number of contacts different in this and refererence complex.
        @rtype: int
        """
        both = N0.logical_or( self.resContacts(cutoff), ref.resContacts(cutoff))
        return N0.sum(N0.sum(both)) - self.contactsShared( ref, cutoff )


    def fractionNativeContacts(self, ref, cutoff=None ):
        """
        Fraction of native B{residue-residue} contacts.
        
        @param ref: native complex
        @type  ref: Complex
        @param cutoff: maximal atom-atom distance, None .. previous setting
        @type  cutoff: float

        @return: fraction of native contacts
        @rtype: float
        """
        cont     = self.resContacts( cutoff, refComplex=ref )
        ref_cont = ref.resContacts( cutoff )

        result = N0.sum(N0.sum( ref_cont * cont ))*1.0
        return result / N0.sum( N0.sum( ref_cont ))


    def fractionNativeSurface(self, cont, contRef ):
        """
        fraction of atoms/residues that are involved in B{any} contacts
        in both complexes.

        @param cont: contact matrix
        @type  cont: matrix
        @param contRef: reference contact matrix
        @type  contRef: matrix
        
        @return: (fractRec, fractLig), fraction of atoms/residues that
                  are involved in any contacts in both complexes
        @rtype: (float, float)
           
        """
        lig, ligRef = N0.clip( N0.sum(cont),0,1),  N0.clip( N0.sum(contRef), 0,1)
        rec    = N0.clip( N0.sum(cont, 1),0,1)
        recRef = N0.clip( N0.sum(contRef, 1), 0,1)

        fLig = N0.sum( N0.logical_and( lig, ligRef )) *1./ N0.sum( ligRef )
        fRec = N0.sum( N0.logical_and( rec, recRef )) *1./ N0.sum( recRef )

        return (fRec, fLig)


    def rmsLig( self, ref ):
        """
        Rms of ligand from reference ligand after superimposing the two
        receptors.
        
        @param ref: reference complex, must have identical atoms
        @type  ref: Complex
        
        @return: ligand rmsd
        @rtype: float

        @note: not implemented
        """
        pass


    def rmsInterface( self, ref, cutoff=4.5, fit=1 ):
        """
        Rmsd between this and reference interface. The interface is
        defined as any residue that has an atom which is within the
        distance given by |cutoff| from its partner.
        
        @param ref: reference complex
        @type  ref: Complex
        @param cutoff: atom distance cutoff for interface residue definition
                       (default: 4.5)
        @type  cutoff: float
        @param fit: least-squares fit before calculating the rms (default: 1)
        @type  fit: 1|0
        
        @return: interface rmad
        @rtype: float
        """
        ## casting
        this = self
        if not ref.rec_model.equals( self.rec_model )[1] \
           or not ref.lig_model.equals( self.lig_model )[1]:

            m_rec, m_rec_ref, m_lig, m_lig_ref = self.equalAtoms( ref )
            this = self.compress( m_rec, m_lig )
            ref  = ref.compress( m_rec_ref, m_lig_ref )

        ## determine interface
        contacts = ref.resContacts( cutoff )

        if_rec = ref.rec_model.res2atomMask( N0.sum( contacts, 1 ) )
        if_lig = ref.lig_model.res2atomMask( N0.sum( contacts, 0 ) )

        mask_interface = N0.concatenate( (if_rec, if_lig) )
        mask_heavy = N0.concatenate( (ref.rec().maskHeavy(),
                                   ref.lig_model.maskHeavy()) )
        mask_interface = mask_interface * mask_heavy

        ## rms
        ref_model = ref.model()
        this_model= this.model()

        return ref_model.rms( this_model, mask_interface, fit=fit)


    def contactResPairs(self, cm=None ):
        """
        Get list of residue names paired up in contacts.
        
        @param cm: pre-calculated contact matrix (default: None)
        @type  cm: matrix
        
        @return: list of tuples [('N','G'), ('P','C')..]
        @rtype: [(str.str)..]
        """
        if cm is None:
            cm = self.resContacts()

        seq_lig = self.lig().sequence()
        seq_rec = self.rec().sequence()

        result = []

        for i in range( 0, len(seq_rec) ):
            recRes = seq_rec[i]

            for j in range( 0, len(seq_lig)):
                ligRes = seq_lig[j]

                if cm[i,j]:
                    result += [(recRes, ligRes)]

        return result


    def contactResDistribution( self, cm=None ):
        """
        Count occurrence of residues in protein-protein interface.
        
        @param cm: pre-calculated contact matrix (default: None)
        @type  cm: matrix
        
        @return: dict {'A':3, 'C':1, .. } (20 standard amino acids)
        @rtype: dict
        """
        if cm is None:
            cm = self.resContacts()

        ## get mask for residues involved in contacts
        maskLig = N0.sum( cm )
        maskRec = N0.sum( N0.transpose( cm ))

        ## get sequence of contact residues only
        seqLig = N0.compress( maskLig, self.lig().sequence() )
        seqRec = N0.compress( maskRec, self.rec().sequence() )
        seq    = ''.join( seqLig ) + ''.join(seqRec) ## convert back to string

        ## count occurrence of letters
        result = {}
        for aa in molUtils.allAA():
            result[aa] = seq.count( aa )

        return result


    def contactTypeDistribution( self, cm = None ):
        """
        Count occurrence of residue types aromatic (a), charged(c),
        polar(p), else/unpolar(u).
        
        @param cm: pre-calculated contact matrix (default: None)
        @type  cm: matrix
        
        @return: dict {'a':4, 'c':5, 'p':12, 'u':4}
        @rtype: dict
        """
        resDistr = self.contactResDistribution( cm )

        result = {'a':0, 'c':0, 'p':0, 'u':0}

        for r in resDistr.keys():
            result[ molUtils.resType( r )[0] ] += resDistr[ r ]

        ## normalize
        s = N0.sum( result.values() )
        for r in result.keys():
            result[r] = result[r]*1.0 / s

        return result


    def resPairCounts(self, cm=None ):
        """
        Count how often (all possible) res-res combinations occur
        in the given contacts.
        
        @param cm: pre-calculated contact matrix  (default: None)
        @type  cm: matrix 
        
        @return: dict {'AA':3,'AC':0, .. 'WW':2, 'WY':0 }
        @rtype: dict
        """
        resPairs = self.contactResPairs( cm )
        allAA = molUtils.allAA()

        allPairs = {}

        ## create dictionary entries for each possible pair of residues
        for i in range(0, len(allAA)):

            for j in range(i, len( allAA ) ):
                key = t.sortString( allAA[j]+allAA[i] )
                allPairs[key] = 0

        ## count occurrence of each pair
        for pair in resPairs:

            key = t.sortString( pair[0] + pair[1] )
            allPairs[ key ] += 1

        return allPairs


    def loadResContacts( self ):
        """
        Uncompress residue contact matrix if necessary.
        
        @return: dict with contact matrix and parameters OR None
        @rtype: dict OR None
        """
        ## Backwards compatibility
        if self.contacts is not None and type( self.contacts ) == str:
            self.contacts = t.load( self.contacts )
            EHandler.warning("loading old-style pickled contacts.") 
            return self.contacts

        ## New, uncompression from list of indices into raveled array
        if self.contacts is not None and \
           len( N0.shape( self.contacts['result'])) == 1:

            try:
                lenRec, lenLig = self.contacts['shape']
            except:
                EHandler.warning("uncompressing contacts without shape")
                lenRec = self.rec().lenResidues()
                lenLig = self.lig().lenResidues()

            m = N0.zeros( lenRec * lenLig )
            N0.put( m, self.contacts['result'], 1 )

            self.contacts['result'] = N0.reshape( m, (lenRec, lenLig) )

        return self.contacts


    def resContacts(self, cutoff=4.5, maskRec=None, maskLig=None,
                    refComplex=None, force=0, cache=1, cache_pw=0 ):
        """
        Matrix of all residue - residue contacts between receptor and
        ligand. Result is cached.
        
        @param cutoff: float/int, cutoff in \\AA for atom-atom contact to be
                       counted ( default 4.5; if None, last one used or 4.5)
        @type  cutoff: float
        @param maskRec: atom mask, receptor atoms to consider
        @type  maskRec: [1|0]
        @param maskLig: atom mask, ligand atoms to consider
        @type  maskLig: [1|0]
        @param force: re-calculate even if cached matrix is available
                      (default: 0)
        @type  force: 0|1
        @param cache: cache result (default: 1)
        @type  cache: 0|1
        @param cache_pw: cache pairwise atom distances (default:0)
        @type  cache_pw: 0|1
        @param refComplex: if <> None, 'reshape' contact matrix, so that
                           residue positions match residue positions in
                           matrices from castComplex. Useful if calling
                           complex is a crystallized reference structure
                           with slight sequence variations from the docked
                           receptor and ligand.
        @type  refComplex: Complex

        @return: residue contact matrix,
                 2-D array(residues_receptor x residues_ligand) of 0 or 1
        @rtype: array
        """
        if cutoff is None:
            if self.contacts is not None:
                cutoff = self.contacts['cutoff']
            else:
                cutoff = 4.5

        if not force and self.loadResContacts() is not None \
           and self.contacts['cutoff'] == cutoff \
           and self.contacts['maskRec'] == maskRec \
           and self.contacts['maskLig'] == maskLig:

            result = self.contacts['result']

        else:
            result = self.__resContacts( cutoff, maskRec, maskLig, cache_pw )

        if cache:
            self.contacts = {}
            self.contacts['cutoff'] = cutoff
            self.contacts['maskRec'] = maskRec
            self.contacts['maskLig'] = maskLig
            self.contacts['result'] = result

        # delete/insert rows or columns to match sequence of reference complex
        if refComplex is not None:
            result = self.__alignMatrixDimension(result,
                                                 self.rec_model.sequence(),
                                                 refComplex.rec_model.sequence() )
            result = self.__alignMatrixDimension(result,
                                                 self.lig_model.sequence(),
                                                 refComplex.lig_model.sequence(),
                                                 1)

        return result


    def __alignMatrixDimension(self, cm, thisSeq, castSeq, axis=0):
        """
        Correct one dimension of contactMatrix by inserting and deleting
        columns, so that it can be later compared to contact matrices based
        on slightly different sequences.
        
        @param cm: contact matrix, 2D matrix of residue contacts
                   recceptor x ligand sequence
        @type  cm: array
        @param thisSeq: AA sequence of this dimension of the contactMatrix
        @type  thisSeq: string
        @param castSeq: AA sequence of this dimension in the other contact
        @type  castSeq: string
        @param axis: which dimension to adapt (0=receptor, 1=ligand)
        @type  axis: 1|0
        
        @return: contact matrix with residue contacts compatible to refSeq.
        @rtype: 2D array
        """
        # compare the two sequences
        seqdiff = SequenceMatcher(None, thisSeq, castSeq)
        seqDiff = seqdiff.get_opcodes()
        ## print seqDiff

        # decide which dimension to work on
        if not axis:
            cm = N0.transpose( cm )

        seqCount = 0   # keep track of sequence length changes
        i=0

        for list in seqDiff:

            # remove the column corresponding to the deletion in the
            # docked sequence
            if str( seqDiff[i][0] ) == 'delete':

                # separate matrix into before and after deletion
                matrixSeg1 = cm[ :, : seqDiff[i][1] + seqCount ]
                matrixSeg2 = cm[ :, seqDiff[i][2] + seqCount : ]
                # concatenate part
                cm = N0.concatenate( ( matrixSeg1, matrixSeg2 ), 1)
                seqCount = seqCount + seqDiff[i][1] - seqDiff[i][2]

            # inserts zeros in the column where there is a insertion in the
            # docked sequence
            if str( seqDiff[i][0] ) == 'insert':

                # create a matrix to be inserted
                insertZeros= seqDiff[i][4] - seqDiff[i][3]
                insertColumns = N0.array( [ [0] * insertZeros ] * N0.size(cm,0) )
                # separate matrix into before and after insertion
                matrixSeg1 = cm[ :, : seqDiff[i][1] + seqCount ]
                matrixSeg2 = cm[ :, seqDiff[i][2] + seqCount : ]
                # concatenate parts with the zero matrix
                cm = N0.concatenate( (matrixSeg1,insertColumns,matrixSeg2), 1)
                seqCount = seqCount + seqDiff[i][4] - seqDiff[i][3]

            i=i+1

        if not axis:
            return N0.transpose( cm )
        return cm


    def __unmaskedMatrix( self, contacts, rec_mask, lig_mask ):
        """
        Map contacts between selected rec and lig atoms back to all atoms
        matrix.
        
        @param contacts: contact matrix, array sum_rec_mask x sum_lig_mask
        @type  contacts: array
        @param rec_mask: atom mask
        @type  rec_mask: [1|0]
        @param lig_mask: atom mask
        @type  lig_mask: [1|0]

        @return: atom contact matrix, array N_atoms_rec x N_atoms_lig
        @rtype: array
        """
        l_rec = len( self.rec_model )
        l_lig = len( self.lig_model )

        ## map contacts back to all atoms matrix
        r = N0.zeros( l_rec * l_lig )
        rMask = N0.ravel( N0.outerproduct( rec_mask, lig_mask ) )

        ## (Optimization: nonzero is time consuming step)
        N0.put( r, N0.nonzero( rMask ), N0.ravel( contacts ) )

        return N0.resize( r, (l_rec, l_lig))


    def __atomContacts(self, cutoff, rec_mask, lig_mask, cache):
        """
        Intermolecular distances below cutoff after applying the two masks.
        
        @param cutoff: cutoff for B{atom-atom} contact in \\AA
        @type  cutoff: float
        @param rec_mask: atom mask
        @type  rec_mask: [1|0]
        @param lig_mask: atom mask
        @type  lig_mask: [1|0]
        @param cache: cache pairwise atom distance matrix
        @type  cache: 1|0
        
        @return: atom contact matrix, array sum_rec_mask x sum_lig_mask
        @rtype: array
        """
        ## get atom coordinats as array 3 x all_atoms
        rec_xyz = self.rec().getXyz()
        lig_xyz = self.lig().getXyz()

        ## get pair-wise distances -> atoms_rec x atoms_lig
        dist = getattr( self, 'pw_dist', None )
        if dist is None or \
               N0.shape( dist ) != ( N0.sum(rec_mask), N0.sum(lig_mask) ):
            dist = self.__pairwiseDistances(N0.compress( rec_mask, rec_xyz, 0),
                                            N0.compress( lig_mask, lig_xyz, 0) )
        if cache:
            self.pw_dist = dist

        ## reduce to 1 (distance < cutoff) or 0 -> n_atoms_rec x n_atoms_lig
        return N0.less( dist, cutoff )


    def atomContacts( self, cutoff=4.5, rec_mask=None, lig_mask=None, cache=0,
                      map_back=1 ):
        """
        Find all inter-molecular B{atom-atom} contacts between rec and lig
        
        @param cutoff: cutoff for atom - atom contact in \\AA
        @type  cutoff: float
        @param rec_mask: atom mask (default: all heavy)
        @type  rec_mask: [1|0]
        @param lig_mask: atom mask (default: all heavy)
        @type  lig_mask: [1|0]
        @param cache: cache pairwise atom distance matrix (default: 0)
        @type  cache: 1|0
        @param map_back: map masked matrix back to matrix for all atoms
        @type  map_back: 1|0
        
        @return: atom contact matrix, Numpy array N(atoms_lig) x N(atoms_rec)
        @rtype: array
        """
        if lig_mask is None:
            lig_mask = self.lig().maskHeavy()

        if rec_mask is None:
            rec_mask = self.rec().maskHeavy()

        contacts = self.__atomContacts( cutoff, rec_mask, lig_mask, cache )

        if not map_back:
            ## contact matrix after masking rec and lig
            return contacts

        return self.__unmaskedMatrix( contacts, rec_mask, lig_mask )


    def __resContacts(self, cutoff, maskRec=None, maskLig=None, cache=0 ):
        """
        Find all inter-molecule B{residue-residue} contacts between receptor
        and ligand. A contact between A and B is set if any heavy atom of
        A is within |cutoff| A of B.
        
        @param cutoff: distance cutoff in \\AA
        @type  cutoff: float
        @param maskRec: atom mask (default: all heavy)
        @type  maskRec: [1|0]
        @param maskLig: atom mask (default: all heavy)
        @type  maskLig: [1|0]
        @param cache: cache pairwise atom distance matrix to pw_dist
                      (default:0)
        @type  cache: 1|0
        
        @return: residue contact matrix, 2-D Numpy
                 N0.array(residues_receptor x residues_ligand) where
                 1-contact, 0-no contact
        @rtype: array
        """
        ## get contact matrix atoms_rec x atoms_lig
        c = self.atomContacts( cutoff, maskRec, maskLig, cache )

        ## convert atoms x atoms to residues x residues matrix
        return self.__atom2residueMatrix( c )


    def __atom2residueMatrix( self, m ):
        """
        Reduce binary matrix of n x k atoms to binary matrix of i x j residues.
        
        @param m: atom contact matrix,
                  array n x k with 1(contact) or 0(no contact)
        @type  m: array
        
        @return: residue contact matrix,
                 2-D numpy array(residues_receptor x residues_ligand)
        @rtype: array
        """
        recInd = N0.concatenate((self.rec().resIndex(),
                              [ self.rec().lenAtoms()] ))
        ligInd = N0.concatenate((self.lig_model.resIndex(),
                              [ self.lig_model.lenAtoms() ] ))

        residueMatrix = N0.zeros(( len(recInd)-1, len(ligInd)-1 ), N0.Int)

        for r in range( len(recInd)-1 ):

            for l in range( len(ligInd)-1 ):

                res2res = m[ int(recInd[r]):int(recInd[r+1]),
                             int(ligInd[l]):int(ligInd[l+1]) ]

                if N0.any( res2res ):
                    residueMatrix[r, l] = 1

        return residueMatrix


    def equalAtoms( self, ref ):
        """
        Compare two complexes' atom content of receptor and ligand. 
        Apply to SORTED models without HETATOMS. Coordinates are not checked.

        @param ref: reference complex
        @type  ref: Complex

        @return: (mask_rec, mask_lig, mask_rec_ref, mask_lig_ref),
                 4 atom masks for all equal atoms
        @rtype: [1|0], [1|0], [1|0], [1|0]
           
        @note: in some rare cases m1.equalAtoms( m2 ) gives a different
               result than m2.equalAtoms( m1 ). This is due to the used
               SequenceMatcher class.           
        """
        m_rec, m_rec_ref = self.rec_model.equalAtoms( ref.rec_model )
        m_lig, m_lig_ref = self.lig_model.equalAtoms( ref.lig_model )
        return m_rec, m_lig, m_rec_ref, m_lig_ref


    def compareAtoms( self, ref ):
        """
        Compare atom content and oder of rec and lig of 2 complexes.
        Get list of all atom positions of rec and lig in both complexes that
        have a matching atom (by residue and atom name) in the other complex.
        Indices for this complex are returned in the right oder to match the
        atom order of ref. I.e., in contrast to equalAtoms, castAtoms can
        equalize two complexes where atoms are ordered differently within the
        residues.

        @param ref: reference complex
        @type  ref: Complex
        
        @return: (ind_rec, ind_lig, ind_rec_ref, ind_lig_ref),
                  4 lists of indices
        @rtype: [int], [int], [int], [int]
        """
        i_rec, i_rec_ref = self.rec_model.compareAtoms( ref.rec_model )
        i_lig, i_lig_ref = self.lig_model.compareAtoms( ref.lig_model )

        return i_rec, i_lig, i_rec_ref, i_lig_ref


    def take( self, rec_pos, lig_pos ):
        """
        Get copy of this complex with given atoms of rec and lig.

        @param rec_pos: receptor indices to take
        @type  rec_pos: [int]
        @param lig_pos: ligand  indices to take
        @type  lig_pos: [int]

        @return: new complex
        @rtype: Complex
        """
        r = self.__class__()
        r.lig_model = self.lig_model.take( lig_pos )
        r.rec_model = self.rec_model.take( rec_pos )
        r.info = deepcopy( self.info )

        if self.pw_dist:
            r.pw_dist = N0.take( self.pw_dist, rec_pos, 1 )
            r.pw_dist = N0.take( r.pw_dist, lig_pos )

        r.ligandMatrix = copy( self.ligandMatrix )

        ## todo: take cached contacts as well

        return r


    def compress( self, rec_mask, lig_mask ):
        """
        Compress complex using a rec and lig mask.
        
        @param rec_mask: atom mask 
        @type  rec_mask: [1|0]
        @param lig_mask: atom mask 
        @type  lig_mask: [1|0]

        @return: compressed complex
        @rtype: Complex
        """
        return self.take( N0.nonzero( rec_mask ), N0.nonzero( lig_mask ) )


    def contPairScore(self, cutoff=6.0):
        """
        Score interaction surface residue pairs.
        Info on Scoring matrix see L{Biskit.molUtils}

        @param cutoff: CB-CB distance cutoff for defining a contact
                       (default: 6.0)
        @type  cutoff: float

        @return: score
        @rtype: float
        """
        score = 0
        cm = self.resContacts(cutoff,self.rec().maskCB(), self.lig().maskCB(),
                              cache=0 )

        pairFreq = self.resPairCounts(cm)

        for pair in pairFreq:
            score += pairFreq[pair]*molUtils.pairScore[pair]

        return score


    def interfaceArea( self, profiles=0, log=StdLog(), verbose=1 ):
        """
        Calculate the difference between the surface area of the
        complex vs. its free components in square angstrom.

        @param profiles: option to return the lig, rec and com profiles
                          rather than the value (both absolute and relative
                          values are returned)
        @type  profiles: 1|0 (default: 0)
        @param log: log file [STDOUT]
        @type  log: Biskit.LogFile
        @param verbose: give progress report [1]
        @type  verbose: bool | int

        @return: AS area, MS area OR
                 a dictionary of lig, rec and com profiles
        @rtype: (float, float) or dict      
        """
        rcom = self.rec()
        lcom = self.lig()
        ccom = self.model()
        result = {}

        def getSurface( model, key ):
            if verbose:
                log.write("Calculating SurfaceRacer data for %s..."%key)
                d = PDBDope( model )
                d.addSurfaceRacer( probe=1.4 )
            if verbose:
                log.writeln('Done.')
                result['%s_AS'%key] = model.profile('AS')
                result['%s_MS'%key] = model.profile('MS')
                result['%s_relAS'%key] = model.profile('relAS')
                result['%s_relMS'%key] = model.profile('relMS')

        getSurface( rcom, 'rec' )
        getSurface( lcom, 'lig' )
        getSurface( ccom, 'com' )

        if not profiles:
            return sum(result['rec_AS']) + sum(result['lig_AS']) - \
                   sum(result['com_AS']),\
                   sum(result['rec_MS']) + sum(result['lig_MS']) - \
                   sum(result['com_MS'])              
        else:
            return result
        

##     def prosaProfile( self ):
##         """
##         Call ProsaII and calculate PROSA energy profile for the Complex.
##         -> array 3 x N_atoms (pair energy, surface energy, combined energy )
##         """
##         m = self.model()

##         ## make PROSA compatible
##         m.renumberResidues()
##         for a in m.getAtoms():
##             a['chain_id'] = 'P'

##         # write temp pdb to disk
##         f = tempfile.mktemp('prosa_pdb')
##         m.writePdb( f, ter=0)

##         try:
##             # run prosa using Wolfgangs Prosa class
##             prosa = ProsaII(executable = prosaBin, temp_dir = t.tempDir() )
##             energies = prosa.analyseEnergy( f )
##         except:
##             print 'Prosa result could not be parsed. CHECK PROSA EXECUTABLE!'

##         # delete temp pdb file 
##         os.system('rm '+ f)

##         return energies


##     def prosaEnergy( self ):
##         """
##         Calculate ProsaII total energies for the complex.
##         -> N0.array(pair energy, surface energy, combined energy ) 
##         """
##         return N0.sum( self.prosaProfile() )


    def prosa2003Energy( self ):
        """
        Calculate Prosa2003 total energies for the complex.
        
        @return: Prosa energy profiles,
                 N0.array(pair energy, surface energy, combined energy )
        @rtype: (array, array, array)
        """
        rec, lig = self.rec_model, self.lig_model
        p = Prosa2003( [rec, lig], debug=1 )
        p.run()

        return p.prosaEnergy()


    def foldXEnergy( self, force=1 ):
        """
        Calculate E_rec and/or E_lig only if not given::
          E_com - (E_rec + E_lig).
        
        @param force: force calc of E_rec and E_lig
        @type  force: 1|0
        
        @return: dict with the different fold-X energy terms
        @rtype: dict
        """
        rec, lig = self.rec_model, self.lig_model

        ## add/update lig/rec fold-X energies if necessary
        if not 'foldX' in rec.info or rec.xyzChanged or force:
            d = PDBDope( rec )
            d.addFoldX()

        if not 'foldX' in lig.info or lig.xyzChanged or force:
            d = PDBDope( lig )
            d.addFoldX()
            try:
                self.lig_transformed.info = lig.info
            except:
                pass
            
        m = self.model()

        d = PDBDope( m )
        d.addFoldX()

        e_rec = rec.info['foldX']
        e_lig = lig.info['foldX']
        e_com = m.info['foldX']

        r = {}
        for key in e_com:
            e = e_com[key] - ( e_lig[key] + e_rec[key] )
            r[key] = e

        return r


    def conservationScore( self, cons_type='cons_ent', ranNr=150,
                           log=StdLog(), verbose=1 ):
        """
        Score of conserved residue pairs in the interaction surface.
        Optionally, normalized by radom surface contacts.

        @param cons_type: precalculated conservation profile name,
                          see L{Biskit.PDBDope}.
        @type  cons_type: str
        @param ranNr: number of random matricies to use (default: 150)
        @type  ranNr: int
        @param log: log file [STDOUT]
        @type  log: Biskit.LogFile
        @param verbose: give progress report [1]
        @type  verbose: bool | int

        @return: conservation score
        @rtype: float
        """
        try:
            recCons = self.rec().profile( cons_type, updateMissing=1 )
        except:
            if verbose:
                log.add('\n'+'*'*30+'\nNO HHM PROFILE FOR RECEPTOR\n'+\
                        '*'*30+'\n')
            recCons = N0.ones( self.rec().lenResidues() )
        try:
            ligCons = self.lig().profile( cons_type, updateMissing=1 )
        except:
            if verbose:
                log.add(\
                            '\n'+'*'*30+'\nNO HHM PROFILE FOR LIGAND\n'+'*'*30+'\n')
            ligCons = N0.ones( self.lig().lenResidues() )

        if self.rec().profile( 'surfMask' ):
            recSurf = self.rec().profile( 'surfMask' )
        else:
            d = PDBDope(self.rec())
            d.addSurfaceMask()

        if self.lig().profile( 'surfMask' ):
            ligSurf = self.lig().profile( 'surfMask' )
        else:
            d = PDBDope(self.lig())
            d.addSurfaceMask()

        surfMask = N0.ravel(N0.outerproduct( recSurf, ligSurf ))

        missing = N0.outerproduct( N0.equal( recCons, 0), N0.equal(ligCons,0))

        cont = self.resContacts() * N0.logical_not(missing)

        consMat = N0.outerproduct( recCons, ligCons )

        score = cont* consMat

        # get a random score
        if ranNr != 0:
            if self.verbose:
                self.log.write('.')
            ranMat =  mathUtils.random2DArray( cont, ranNr, mask=surfMask )
            random_score = N0.sum(N0.sum( ranMat * consMat ))/( ranNr*1.0 )
            return N0.sum(N0.sum(score))/random_score

        else:
            return N0.sum(N0.sum(score))/ N0.sum(N0.sum(cont))


    def rtTuple2matrix( self, r, t):
        """
        Put rotation and translation matrix into single 4x4 matrix.
        
        @param r: rotation matric, array 3x3 of float
        @type  r: array
        @param t: translation vector, array 1x3 of float
        @type  t: vector
        
        @return: rotation/translation matrix, array 4x4 of float
        @rtype: array
        """
        ## create 3 x 4 matrix: 0:3, 0:3 contains rot; 3,0:3 contains trans
        result = N0.concatenate( (r, N0.transpose( [ t.tolist() ] )), 1)
        ## make it square
        result = N0.concatenate( (result, N0.array([[0,0,0,1]],N0.Float32)), 0)

        return result.astype(N0.Float32)


    def extractLigandMatrix( self, lig):
        """
        Find transformation matrix for rigid body-transformed ligand.

        @param lig: ligand model
        @type  lig: PDBModel
        
        @return: rotation/translation matrix, array 4x4 of float
        @rtype: array
        """
        lig = lig.compress( lig.maskCA() )
        template = self.lig_model.compress( self.lig_model.maskCA() )

        r,t = self.__findTransformation( lig.xyz, template.xyz )

        return self.rtTuple2matrix( r, t )


    def __findTransformation(self, x, y):
        """
        Match two arrays by rotation and translation. Returns the
        rotation matrix and the translation vector.
        Back transformation:
        for atom i new coordinates will be::
            y_new[i] = N0.dot(r, y[i]) + t
            
        for all atoms in one step::
            y_new = N0.dot(y, N0.transpose(r)) + t

        @param x: coordinates
        @type  x: array
        @param y: coordinates
        @type  y: array

        @return: rotation matrix, translation vector
        @rtype: array, array      
        
        @author: Michael Habeck
        """
        from numpy.linalg import svd

        ## center configurations
        x_av = N0.sum(x) / len(x)
        y_av = N0.sum(y) / len(y)
        x = x - x_av
        y = y - y_av
        ## svd of correlation matrix
        v, l, u = svd(N0.dot(N0.transpose(x), y))
        ## build rotation matrix and translation vector
        r = N0.dot(v, u)
        t = x_av - N0.dot(r, y_av)

        return r, t


    def __pairwiseDistances(self, u, v):
        """
        pairwise distance between 2 3-D numpy arrays of atom coordinates.

        @param u: coordinates
        @type  u: array
        @param v: coordinates
        @type  v: array
        
        @return: Numpy array len(u) x len(v)
        @rtype:array
        
        @author: Wolfgang Rieping.
        """
        ## check input
        if not type( u ) == arraytype or\
           not type( v ) == arraytype:
            raise ComplexError('unsupported argument type ' + \
                               str( type(u) ) + ' or ' + str( type(v) ) )

        diag1= N0.diagonal(N0.dot(u,N0.transpose(u)))
        diag2= N0.diagonal(N0.dot(v,N0.transpose(v)))
        dist= -N0.dot(v,N0.transpose(u))-N0.transpose(N0.dot(u,N0.transpose(v)))
        dist= N0.transpose(N0.asarray(list(map(lambda column,a:column+a, \
                                   N0.transpose(dist), diag1))))

        return N0.transpose(N0.sqrt(N0.asarray(
            list(map(lambda row,a: row+a, dist, diag2)))))


    ## obsolete
    def __extractLigandMatrix(self, fcomplex):
        """
        Compare structure from hex complex with original ligand pdb
        and store transformation matrix of ligand in self.ligandMatrix.
        
        @param fcomplex: pdb file with hex complex
        @type  fcomplex: complec
        
        @return: rotation matrix and translation matrix as tuple
        @rtype: (array, array)
        """
        docked_pdb = self._extractLigandStructure(fcomplex)

        xyz_docked = N0.compress( docked_pdb.maskCA(), docked_pdb.xyz )
        xyz_template = N0.compress( self.lig_model.maskCA(),
                                 self.lig_model.xyz )

        (r, t) = self._findTransformation(xyz_docked, xyz_template)
        return (r,t)



#############
##  TESTING        
#############
import biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test case"""

    def test_Complex(self):
        """Dock.Complex test"""
        from biskit import XplorModel

        lig = XplorModel( t.testRoot() + "/com/1BGS.psf",
                        t.testRoot() + "/com/lig.model")

        rec = XplorModel( t.testRoot() + "/com/1BGS.psf",
                        t.testRoot() + "/com/rec.model")

        rec = rec.compress( rec.maskHeavy() )
        lig = lig.compress( lig.maskHeavy() )

        c = Complex(rec, lig)
        c.info['soln'] = 1

        cont = c.atomContacts( 6.0 )
        contProfile_lig = N0.sum( cont )
        contProfile_rec = N0.sum( cont, 1 )

        try:
            dope = PDBDope( c.rec_model )
            dope.addSurfaceRacer( probe=1.4 )
            rec_surf = c.rec_model.profile2mask( 'MS', 0.0000001, 1000 )

            dope = PDBDope( c.lig_model )
            dope.addSurfaceRacer( probe=1.4 )
            lig_surf = c.lig_model.profile2mask( 'MS', 0.0000001, 1000 )
        except:
            pass

        if self.local:           
            from biskit.exe import Pymoler

            self.pm = Pymoler()
            self.pm.addPdb( c.rec(), 'rec' )
            self.pm.addPdb( c.lig(), 'lig' )

            self.pm.colorAtoms( 'rec', contProfile_rec )
            self.pm.colorAtoms( 'lig', contProfile_lig )

            rec_sphere = c.rec().clone()
            rec_sphere.xyz = mathUtils.projectOnSphere( rec_sphere.xyz )

            lig_sphere = c.lig().clone()
            lig_sphere.xyz = mathUtils.projectOnSphere( lig_sphere.xyz )

            self.pm.addPdb( rec_sphere, 'rec_sphere' )
            self.pm.addPdb( lig_sphere, 'lig_sphere' )

            self.pm.colorAtoms( 'rec_sphere', contProfile_rec )
            self.pm.colorAtoms( 'lig_sphere', contProfile_lig )

            self.pm.add( 'hide all')

            self.pm.add( 'color grey, (b=0)' )
            self.pm.add( 'show stick, (rec or lig)' )
            self.pm.add( 'show surf, rec_sphere')

            self.pm.add( 'zoom all' )

            self.pm.show()
            
            globals().update( locals() )

        self.assertEqual( N0.sum(contProfile_lig) + N0.sum(contProfile_rec),
                          2462 )
   

if __name__ == '__main__':

    BT.localTest()

