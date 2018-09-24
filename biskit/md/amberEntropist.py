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
Run ptraj entropy analysis on Trajectory instance.
"""

##import numpy.oldnumeric as oldN
import biskit.core.oldnumeric as N0

import tempfile, os
import random, time

import biskit.tools as t
import biskit.mathUtils as MU
## import Biskit.settings as settings
from biskit.errors import BiskitError
from biskit import PDBModel, EHandler, LocalPath
from biskit.dock import Complex

## allow relative imports when calling module by itself for testing (pep-0366)
if __name__ == "__main__" and __package__ is None:
    import biskit.md; __package__ = "biskit.md"

from .trajectory import Trajectory
from .ensembleTraj import EnsembleTraj
from .amberCrdEntropist import AmberCrdEntropist, EntropistError
from .amberParmBuilder import AmberParmBuilder


class AmberEntropist( AmberCrdEntropist ):
    """
    Run ptraj entropy analysis on Trajectory instance.
    """

    def __init__( self, traj=None, parm=None, crd=None, ref=None, cast=0,
                  chains=None, border=None, split=0, shift=0, shuffle=0,
                  thin=None,
                  s=0, e=None, ss=0, se=None, step=1, atoms=None, heavy=0,
                  solvent=0, protein=0,
                  ex=[], ex_n=0, ex3=None, ex1=None,
                  fit_s=None, fit_e=None, memsave=1,
                  **kw ):
        """
        :param traj: path to 1 or 2 pickled Trajectory instances
                     (2 separated by '+', e.g. 'rec.traj+lig.traj')
        :type  traj: str
        :param parm: try using existing parm file & keep it [create+discard]
        :type  parm: str
        :param crd: target file for amber crd & keep it (default: discard)
        :type  crd: str
        :param ref: superimpose onto this structure
        :type  ref: str|PDBModel|Complex
        :param cast: equalize atom content against ref (if given) (default: 1)
        :type  cast: 0|1

        :param chains: extract chains from traj (default: None, all chains)
        :type  chains: [int]
        :param border: 1st chain of 2nd molecule; required for split, shift,
                       shuffle if traj is not already a tuple of trajectories
        :type  border: int
        :param split: split trajectory after *border* and fit the two halfs
                      separately (default: 0)
        :type  split: 1|0
        :param shift: recombine rec and lig member trajectories, should
                      disrupt correlations between rec and lig, requires
                      *chains* or 2 traj files to identify rec (default: 0)
        :type  shift: int
        :param shuffle: shuffle the order of frames for one trajectory half,
                        requires *border* or 2 traj files to identify rec
        :type  shuffle: 0|1
        :param s: start frame of complete traj (default: 0)
        :type  s: int
        :param e: stop frame of complete traj (default: None)
        :type  e: int
        :param ss: start frame of single member trajectories (only works
                   with EnsembleTraj; overrides s,e) (default: 0)
        :type  ss: int
        :param se: stop frame of single member trajectories (only works
                   with EnsembleTraj; overrides s,e) (default: None)
        :type  se: int                 
        :param step: frame offset (default: 1, no offset)
        :type  step: int
        :param thin: use only randomly distributed fraction of frames
                     (default: all)
        :type  thin: float
        :param atoms: atom names to consider (default: all)
        :type  atoms: [str]
        :param heavy: remove hydrogens (default: 0)
        :type  heavy: 1|0
        :param protein: remove all non-protein atoms (default: don't)
        :type  protein: 1|0
        :param solvent: retain solvent and ions (default: 0)
        :type  solvent: 1|0
        :param ex: exclude member trajectories
        :type  ex: [int] OR ([int],[int])
        :param ex_n: exclude last n members  OR...                
        :type  ex_n: int
        :param ex3: exclude *ex3*rd tripple of trajectories  (default: 0)
                    (index starts with 1! 0 to exclude nothing) OR....
        :type  ex3: int
        :param ex1: exclude *ex1*-th member remaining after applying *ex*
                    (default: None)(index starts with 1! 0 to exclude nothing)
        :type  ex1: int
        :param fit_s: fit to average of different frame slice 
        :type  fit_s: int|None
        :param fit_e: fit to average of different frame slice 
        :type  fit_e: int|None
        :param memsave: delete internal trajectory after writing crd
                        (default: 1)
        :type  memsave: 1|0

        :param kw: additional key=value parameters for AmberCrdEntropist
                   and Executor:
        :type  kw: key=value pairs
        ::
          ... parameters for AmberCrdEntropist
          f_template - str, alternative ptraj input template

          ... and key=value parameters for Executor:
          f_out    - str, target name for ptraj output file (default: discard)
          debug    - 0|1, keep all temporary files (default: 0)
          verbose  - 0|1, print progress messages to log (log != STDOUT)
          node     - str, host for calculation (None->local) NOT TESTED
                          (default: None)
          nice     - int, nice level (default: 0)
          log      - biskit.LogFile, program log (None->STOUT) (default: None)
        """
##         tempfile.tempdir = '/work'

        f_crd = crd or tempfile.mktemp('.crd')
        self.keep_crd = crd is not None

        f_parm = t.absfile( parm or tempfile.mktemp( '.parm' ) )
        self.keep_parm = parm is not None

        self.parmcrd = tempfile.mktemp('_ref.crd')

        AmberCrdEntropist.__init__( self, f_parm, f_crd,
                                    s=s, e=e, step=step, **kw )

        self.fit_s = fit_s ## if None, fit() will set it to 0
        self.fit_e = fit_e ## if None, fit() will set it to len(traj)

        self.cast   = cast
        self.chains = chains
        self.border = border
        self.split  = split
        self.shift  = shift
        self.shuffle= shuffle
        self.thin   = thin
        self.thin_i = None   ## extract same frames from all trajectories
        self.exclude= ex
        self.ex_n   = ex_n

        self.sstart = ss     ## self.start is assigned by AmberCrdEntropist
        self.sstop  = se     ## self.stop  is assigned by AmberCrdEntropist
        self.heavy  = heavy
        self.solvent = solvent
        self.protein = protein
        self.atoms  = atoms
        self.memsave= memsave

        if ex3 in [0, None]:
            self.ex3 = None
        else:
            self.ex3 = ex3-1
        if ex3 is not None: self.ex_n = 0

        if ex1 is not None:
            self.ex3 = None
            self.ex_n = 0
        if ex1 in [0, None]:
            self.ex1 = None
        else:
            self.ex1 = ex1-1

        ## filter atoms
        self.ref = self.prepareRef( ref )

        ## filter atoms, cast to ref, fit to average and ref
        self.traj = self.prepareTraj( traj, self.ref, cast=self.cast )
        self.nframes = len( self.traj )

        ## reset frame limits for AmberCrdEntropist
        self.start = 0
        self.stop = len( self.traj )
        self.step = 1
        

    def prepare( self ):
        """
        Overrides Executor method.
        """
        ## create Parm file for selected atoms of trajectory
        if not os.path.exists( self.f_parm ):
            self.buildParm()
        else:
            if self.verbose: self.log.add('using existing %s' % self.f_parm)

        ## create Amber Crd file for ptraj
        if not os.path.exists( self.f_crd ):
            if self.verbose: self.log.write('Writing amber crd file...')
            self.traj.writeCrd( self.f_crd )
            if self.verbose: self.log.add('done')
            ## release memory
            if self.memsave: self.traj = None
        else:
            if self.verbose: self.log.add('using existing %s' % self.f_crd)


    def buildParm( self ):
        """
        Build amber topology.
        """
        a = AmberParmBuilder( self.traj.ref, verbose=self.verbose,
                              debug=self.debug )
        if self.verbose: self.log.write('Building amber topology...')

        a.parmMirror( self.f_parm, self.parmcrd )
        
        if self.verbose: self.log.add('Topology built')


    def fit( self, traj, refModel=None, mask=None, conv=1e-6 ):
        """
        Fit trajectory until convergence onto it's own average and then
        transform the average of all frames onto the reference.
        
        :param traj: trajectory in which to fit frames
        :type  traj: Trajectory  
        :param refModel: reference PDBModel
        :type  refModel: PDBModel
        :param mask: atom mask for superposition (default: all)
        :type  mask: [1|0]
        :param conv: convergence criteria (default: 1e-6)
        :type  conv: float
        """
        self.fit_e = self.fit_e or len( traj )
        self.fit_s = self.fit_s or 0

        traj.fit( ref=traj.ref, mask=mask, verbose=self.verbose )

        m_avg = traj[self.fit_s : self.fit_e ].avgModel()

        ## fit on average until it's not getting better
        d = 1.
        dd= 1.
        while dd >= conv:

            traj.fit( ref=m_avg, mask=mask, verbose=self.verbose )
            m_new_avg = traj[self.fit_s : self.fit_e].avgModel()

            oldD, d    = d, m_avg.rms( m_new_avg, mask=mask )

            if self.verbose:
                self.log.add( "rms difference: %f" % d )

            dd = oldD - d
            m_avg = m_new_avg

        ## transform trajectory en block onto reference
        if refModel:
            if self.verbose:
                self.log.add('fitting trajectory en-block onto reference...')

            if refModel.atomNames() != traj.ref.atomNames():
                if self.verbose: self.log.add('casting ref for fitting...')
                ref_i, i = refModel.compareAtoms( m_avg )

                refModel = refModel.take( ref_i )
                m_avg    = m_avg.take( i )

                if not mask is None:   mask = N0.take( mask, i )

            r, t = m_avg.transformation( refModel, mask )
            traj.transform( r, t )


    def __splitFilenames( self, f ):
        """
        Split file name::
          split(traj1.dat+traj2.dat) -> (traj1.dat, traj2.dat)
        
        :param f: file name
        :type  f: str

        :return: split filename
        :rtype: str, str
        """
        if f.find("+") != -1 :
            split = f.find("+")
            f1 = f[:split]
            f2 = f[split+1:]

            return f1, f2
        return None


    def __cleanAtoms( self, m ):
        """
        Remove non protein atoms and H if needed.

        :param m: model to clean
        :type  m: PDBModel

        :return: cleaned model
        :rtype: PDBModel      
        """
        if self.protein:            
            m.keep( N0.nonzero( m.maskProtein() ) )
        if self.heavy:
            m.keep( N0.nonzero( m.maskHeavy() ) )
        return m


    def __getModel( self, f ):
        """
        Load PDBModel directly or extract it from Trajectory.
        
        :param f: file name of PDB file, pickled PDBModel, or Trajectory
        :type  f: str
        
        :return: model
        :rtype: PDBModel
        
        :raise IOError: if file does not exist
        """
        p = LocalPath( f )

        isPdb = (f[-4:].upper() == '.PDB' or f[-7:].upper() == '.PDB.GZ')

        if isPdb and p.exists():
            return p.local()

        o = p.load()   ## !! raises IOError if file does not exist

        if isinstance( o, PDBModel ):
            return o
        if isinstance( o, Trajectory ):
            return o.ref

        raise EntropistError('unknown reference type')


    def prepareRef( self, fname ):
        """
        Prepare reference model.
        
        :param fname: file name 
        :type  fname: str

        :return: reference structure
        :rtype: PDBModel|Complex        

        :raise EntropistError: if unknown reference type
        """
        if not fname:
            return None

        if self.__splitFilenames( fname ):
            f1, f2 = self.__splitFilenames( fname )
            m1, m2 = PDBModel( self.__getModel(f1) ), \
                     PDBModel( self.__getModel(f2) )

            ref = Complex( m1, m2 )
        else:
            ref = t.load( fname )

        if isinstance( ref, Trajectory ):
            ref = ref.ref

        if isinstance( ref, PDBModel ):
            return self.__cleanAtoms( ref )

        if isinstance( ref, Complex ):
            self.__cleanAtoms( ref.rec_model )
            self.__cleanAtoms( ref.lig_model )
            ref.lig_model_transformed = None
            return ref

        raise EntropistError('unknown reference type')


    def __add3( self, n_members, excluded, trippleIndex ):
        """
        Add a tripple of numbers from range( n_members ) to be
        excluded for error estimation. Tripples are chosen to have
        minimal overlap. For 10 trajectories (*n_members*=10), the
        first 3 tripples will be (1,2,3), (4,5,6), (7,8,9).

        :param n_members: number of member trajectories
        :type  n_members: int
        :param excluded: excluded member trajectories
        :type  excluded: [ int ]
        :param trippleIndex: 
        :type  trippleIndex: int

        :return: the indices of all excluded member trajectories
        :rtype: [ int ]
        """
        remaining = MU.difference( range( n_members ), excluded )
        tripple = self.tripples( remaining, trippleIndex+1 )[-1]
        return MU.union( excluded, list(tripple) )


    def __add1( self, n_members, excluded, index ):
        """
        Add one number from range( n_members ) to list of excluded indices

        :param n_members: number of member trajectories
        :type  n_members: int
        :param excluded: excluded member trajectories
        :type  excluded: [ int ]
        :param index: 
        :type  index: int

        :return: the indices of all excluded member trajectories
        :rtype: [ int ]
        """
        remaining = MU.difference( range( n_members ), excluded )
        new_i = remaining[index]
        return excluded + [ new_i ]


    def __exclude( self, traj, exclude ):
        """
        Exclude members from a (set of) Trajectory.
        :param traj: input trajectory
        :type traj: EnsembleTraj
        :param exclude: set of indices to be excluded
        :type exclude: [ int ]

        :return: 
        :rtype: EnsembleTraj
        """
        if exclude is None or len( exclude ) == 0:
            return traj

        members = range( traj.n_members )

        if self.verbose:
            self.log.add("excluding members: " + str(exclude))

        return traj.takeMembers( MU.difference( members, exclude ) )


    def __removeMembers( self, t ):
        """
        Some individual trajectories may have to be excluded as
        outliers or for error estimation (depending on the parameters
        passed to AmberEntropist). 

        :param t: one or two ensembles of trajectories
        :type t: EnsembleTraj OR (EnsembleTraj, EnsembleTraj )

        :return: t with some member trajectories excluded, if needed
        :rtype: EnsembleTraj OR (EnsembleTraj, EnsembleTraj )
        """
        if self.ex_n:
            self.exclude = list(range( self.ex_n))

        if type( t ) is tuple and not type( self.exclude ) is tuple:
            self.exclude = ( self.exclude, self.exclude )

        if type( t ) is tuple:
            n_memb = ( t[0].n_members, t[1].n_members )
        else:
            n_memb = t.n_members

        if self.ex3 is not None and type( self.exclude ) is tuple:
            self.exclude = self.__add3(n_memb[0], self.exclude[0], self.ex3),\
                           self.__add3(n_memb[1], self.exclude[1], self.ex3)

        if self.ex3 is not None and type( self.exclude ) is list:
            self.exclude = self.__add3( n_memb, self.exclude, self.ex3 )

        if self.ex1 is not None and type( self.exclude ) is tuple:
            self.exclude = self.__add1(n_memb[0], self.exclude[0], self.ex1),\
                           self.__add1(n_memb[1], self.exclude[1], self.ex1)

        if self.ex1 is not None and type( self.exclude ) is list:
            self.exclude = self.__add1( n_memb, self.exclude, self.ex1 )

        if type( t ) is tuple:
            if not type( self.exclude ) is tuple:
                self.exclude = ( self.exclude, self.exclude )
            t = self.__exclude( t[0], self.exclude[0] ),\
                self.__exclude( t[1], self.exclude[1] )
        else:
            t = self.__exclude( t, self.exclude )

        return t


    def prepareTraj( self, fname, ref=None, cast=1 ):
        """
        Prepare trajectory for Amber.

        :param fname: path to EnsembleTraj OR ( EnsembleTraj, EnsembleTraj )
        :type  fname: str OR (str,str)
        :param ref: reference structure
        :type  ref: EnsembleTraj
        :param cast: cast to reference (same atom content) (default: 1)
        :type  cast: 1|0

        :return: split, fitted or shuffled, etc. trajectory instance
        :rtype: EnsembleTraj OR (EnsembleTraj, EnsembleTraj )
        """
        ## Load 1 or 2
        if self.__splitFilenames( fname ):
            f1, f2 = self.__splitFilenames( fname )
            t = self.loadTraj( f1 ), self.loadTraj( f2 )
        else:
            t = self.loadTraj( fname )
            if self.chains:
                t = t.takeChains( self.chains )

        ## split 1 into 2 if necessary
        if not type(t) is tuple and self.border:
            lig = range( self.border, t.ref.lenChains() )
            t = t.takeChains( range(self.border) ), t.takeChains( lig )

        ## check 2 trajectories were suplied
        if not type(t) is tuple and (self.shift or self.shuffle or self.split):
            raise EntropistError('split,shift,shuffle require -border.')

        ## adapt reference to type of trajectory input 
        if ref and type(t) is tuple and not isinstance( ref, Complex ):
            rec = ref.takeChains( range(t[0].ref.lenChains()) )
            lig = ref.takeChains( range(t[0].ref.lenChains(),
                                        t[1].ref.lenChains()) )
            ref = Complex( rec, lig )

        if ref and type(t) is not tuple and isinstance( ref, Complex ):
            ref = ref.rec_model.concat( ref.lig() )

        ## remove member trajectories (if requested)
        t = self.__removeMembers( t )

        ## cast 1 or 2
        if cast and type( t ) is not tuple:
            self.castTraj( t, ref )

        if cast and type( t ) is tuple:
            self.castTraj( t[0], ref.rec_model )
            self.castTraj( t[1], ref.lig_model )

        ## reorder one half (requires -border or file name pair )
        if self.shift:
            t = self.shiftTraj( t[0], self.shift ), t[1]

        if self.shuffle:
            t = self.shuffleTraj( t[0] ), self.shuffleTraj( t[1] )

        ## fit seperately (requires -border or file name pair)
        if self.split and ref:
            self.fit( t[0], ref.rec() )
            self.fit( t[1], ref.lig() )
        if self.split and not ref:
            self.fit( t[0] )
            self.fit( t[1] )

        if type( t ) is tuple:
            t = t[0].concatAtoms( t[1] )
            ref = ref.rec_model.concat( ref.lig() )

        ## joint fit
        if not self.split:
            self.fit( t, ref )
            
        if self.verbose:
            self.log.add( 'Analysing trajectory with %i atoms and %i frames.' \
                          % (t.lenAtoms(), t.lenFrames()))

        return t


    def load_locked( self, fname ):
        """
        wait with unpickling until another Entropist has finished.

        :param fname: file name
        :type  fname: str

        :return: trajectroy
        :rtype: Trajectroy
        """
        flock = fname + '__locked'

        while os.path.exists( flock ):
            if self.verbose: self.log.write('~')
            time.sleep( random.random() * 10 )
        if self.verbose: self.log.add('')

        try:
            f = open( flock, 'w' )
            f.write('1')
            f.close()

            r = t.load(fname)

        finally:
            t.tryRemove( flock )

        return r


    def loadTraj( self, fname, shift=0 ):
        """
        Load single trajectory.
        """
        if self.verbose:
            self.log.write( 'Loading %s...' % fname)

        traj = self.load_locked( fname )

        if self.verbose:
            self.log.add( 'Processing trajectory...')

        ## convert single member frame index into all member frame index
        if (self.sstart or self.sstop) and isinstance(traj, EnsembleTraj):
            self.start = (self.sstart or 0) * traj.n_members
            self.stop  = (self.sstop  or 0) * traj.n_members
        if (self.sstart or self.sstop) and not isinstance(traj, EnsembleTraj):
            self.start, self.stop = self.sstart, self.sstop
            if self.verbose: self.log.add('Warning: I am using -ss -se instead of -s -e')
        ## remove unwanted frames
        if self.start or self.stop:
            start, stop = self.start, self.stop or len(traj)
            traj = traj[ start : stop ]

        ## stepping (offset)    
        if self.step > 1:
            traj = traj.thin( self.step )

        ## thin with random stepping, use same frames from all trajectories
        if self.thin:
            targetLength = int( round( len( traj ) * self.thin ) )
            self.thin_i = self.thin_i or \
                          MU.randomRange(0, len( traj ), targetLength )
            traj = traj.takeFrames( self.thin_i )

            if self.verbose:
                self.log.add( "Thinned to %i frames." % len( traj ) )

        ## keep only allowed atoms (default: all)
        if self.atoms:
            traj.ref.addChainId()
            aMask = traj.ref.mask( lambda a,ok=self.atoms: a['name'] in ok )
            traj.removeAtoms( N0.nonzero( N0.logical_not( aMask ) )  )

        ## get rid of non-standard atoms, water, ions, etc.
        if not self.solvent:
            l = traj.lenAtoms()
            traj = traj.compressAtoms( N0.logical_not(traj.ref.maskSolvent()) )

            if self.verbose:
                self.log.add('%i solvent/ion atoms deleted.'% (l- traj.lenAtoms()))

        ## delete hydrogens, if requested
        if self.heavy:
            l = traj.lenAtoms()
            traj = traj.compressAtoms( traj.ref.maskHeavy() )
            
            if self.verbose:
                self.log.add('%i hydrogens deleted.' % (l - traj.lenAtoms()) )

        return traj


    def shiftTraj( self, traj, shift=0 ):
        """
        reorder member trajectories
        """
        if not shift:
            return traj

        if not isinstance( traj, EnsembleTraj):
            raise EntropistError('shift requires EnsembleTraj')

        r = list(range( shift, traj.n_members)) + list(range( shift))
        if self.verbose:
            self.log.add('reorder member trajectories: %s...' % str( r ) )

        return traj.takeMembers( r )


    def shuffleTraj( self, traj ):
        """
        reorder all frames at random
        """
        r = list(range( len( traj )))
        random.shuffle( r )
        return traj.takeFrames( r )


    def castTraj( self, traj, refModel ):
        """
        Equalize atom content of traj to refModel. 
        """
        l    = traj.lenAtoms()
        lref = len( refModel )

        if self.verbose:
            self.log.add('comparing traj with %i atoms to reference of %i atoms.'%\
                         (l, lref))
        i, iRef = traj.ref.compareAtoms( refModel )

        refModel.keep( iRef )
        if self.verbose:
            self.log.add("%i atoms deleted from reference."%(lref-len(refModel)))

        traj.keepAtoms( i )
        if self.verbose:
            self.log.add( "%i atoms deleted from trajectory."% (l-len(i) ) )


    def tripples( self, lst, n ):
        """
        Group items of lst into n tripples with minimal overlap.
        """
        all = []
        l = len( lst )

        ## get all possible tripples
        for i in range( l ):
            for j in range( i+1, l ):
                for k in range( j+1, l ):
                    all += [ ( lst[i], lst[j], lst[k] ) ]

        ## calculate pairwise "distance" between tripples
        pw = N0.zeros( (len(all), len(all)), N0.Float32 )
        for i in range( len( all ) ):
            for j in range( i, len(all) ):
                pw[i,j] = pw[j,i] = len( MU.intersection(all[i],all[j]) )**2

        pos = 0
        r = []

        while len( r ) < n:

            r += [ pos ]
            ## overlap of selected tripples with all others
            overlap = N0.sum( N0.array( [ pw[ i ] for i in r ] ) )
            ## select one with lowest overlap to all tripples selected before
            pos = N0.argmin( overlap )

        return N0.take( all, r )


    def cleanup( self ):
        """
        Remove temporary files.
        """
        AmberCrdEntropist.cleanup( self )

        if not self.debug:
            t.tryRemove( self.parmcrd )

            if not self.keep_crd:
                t.tryRemove( self.f_crd )
            if not self.keep_parm:
                t.tryRemove( self.f_parm )


    def finish( self ):
        """
        Called when done.
        """
        AmberCrdEntropist.finish( self )

        ## consistency check of result
        if self.result['nframes'] != self.nframes:
            raise EntropistError('incorrect number of frames: %i instead %i'%\
                  ( self.result['nframes'], self.nframes ))
        
#############
##  TESTING        
#############
import biskit.test as BT

class Test(BT.BiskitTest):
    """Test class"""

    TAGS = [ BT.EXE ]
    
    def test_amberEntropist( self ):
        """AmberEntropist test"""
        import biskit.tools as T
        self.a = AmberEntropist( T.testRoot('/amber/entropy/com_fake.etraj'),
                                 verbose=self.local, debug=self.DEBUG,
                                 log=self.log)
        self.r = self.a.run()
        self.assertTrue( abs(int(self.r['S_total']) - 398) < 2 )
        self.assertAlmostEqual( self.r['mass'], 3254, 0 )
        self.assertAlmostEqual( self.r['S_vibes'], 298, 0 )
        self.assertEqual( int(self.r['S_rot']), 50 )
        self.assertEqual( int(self.r['nframes']), 44 )

if __name__ == '__main__':

    BT.localTest(debug=False)
    
