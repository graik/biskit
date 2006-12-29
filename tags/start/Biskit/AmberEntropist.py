##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##

## last $Author$
## last $Date$
## $Revision$

import tempfile, os, copy
import Numeric as N
import random, time

import Biskit.tools as t
import Biskit.mathUtils as MU
import Biskit.settings as settings
from Biskit.Errors import BiskitError
from Biskit.AmberCrdEntropist import AmberCrdEntropist, EntropistError
from Biskit.AmberParmBuilder import AmberParmBuilder
from Biskit.PDBModel import PDBModel
from Biskit.Trajectory import Trajectory
from Biskit.EnsembleTraj import EnsembleTraj
from Biskit.LocalPath import LocalPath
from Biskit.Dock.Complex import Complex
from Biskit import EHandler

class AmberEntropist( AmberCrdEntropist ):
    """
    Run ptraj entropy analysis on Trajectory instance.
    """

    def __init__( self, traj=None, parm=None, crd=None, ref=None, cast=0,
                  chains=None, border=None, split=0, shift=0, shuffle=0,
                  thin=None,
                  s=0, e=None, ss=0, se=None, step=1, atoms=None, heavy=0, 
                  ex=[], ex_n=0, ex3=None, ex1=None,
                  fit_s=None, fit_e=None, memsave=1,
                  **kw ):
        """
        traj  - str, path to 1 or 2 pickled Trajectory (2 seperated by '+')
        parm  - str, try using existing parm file & keep it [create+discard]
        crd   - str, target file for amber crd & keep it           [discard]
        ref   - str|PDBModel|Complex, superimpose onto this structure
        cast  - 0|1, equalize atom content against ref (if given)       [1]

        chains  - [ int ], extract chains from traj        [None..all chains]
        border  - int, 1st chain of 2nd molecule - required for split, shift,
                  shuffle if traj is not already a tuple of trajectories [] 
        split   - 1|0, split trajectory after |border| and fit the two halfs
                  separately                                   [0]
        shift   - int, recombine rec and lig member trajectories, should
                  disrupt correlations between rec and lig,
                  requires |chains| or 2 traj files to identify rec     [0]
        shuffle - 0|1, shuffle the order of frames for one trajectory half,
                  requires |border| or 2 traj files to identify rec
        s,e     - int, start and stop frame of complete traj.   [0, to end]
        ss, se  - int, start and stop frame of single member trajectories
                  (only works with EnsembleTraj; overrides s,e)
        step    - int, frame offset                                 [no offset]
        thin    - float, use only randomly distributed fraction of frames [all]
        atoms   - [ str ], atom names to consider                       [all]
        heavy   - 1|0, remove hydrogens                                 [0]
        ex      - [int] OR ([int],[int]), exclude member trajectories   [[]]
        ex_n    - int, exclude last n members  OR...                    []
        ex3     - int, exclude |ex3|rd tripple of trajectories          [0]
                  (index starts with 1! 0 to exclude nothing) OR....
        ex1     - int, exclude ex1-th member remaining after applying ex [None]
                  (index starts with 1! 0 to exclude nothing)
        fit_s/e - int | None, fit to average of different frame slice   [=s/e]
        memsave - 1|0, delete internal trajectory after writing crd     [1]

        ... parameters for AmberCrdEntropist
        f_template - str, alternative ptraj input template [default]

        ... and key=value parameters for Executor:
        f_out    - str, target name for ptraj output file        [discard]
        node     - str, host for calculation (None->local)              []
        nice     - int, nice level                                     [0]
        log      - Biskit.LogFile, program log (None->STOUT)            []
        debug    - 0|1, keep all temporary files                       [0]
        verbose  - 0|1, print progress messages to log     [log != STDOUT]
        """
        tempfile.tempdir = '/work'
        
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
        self.atoms  = atoms
        self.memsave= memsave

        if ex3 in [0, None]:
            self.ex3 = None
        else:
            self.ex3 = ex3-1
        if ex3 != None: self.ex_n = 0

        if ex1 != None:
            self.ex3 = None
            self.ex_n = 0
        if ex1 in [0, None]:
            self.ex1 = None
        else:
            self.ex1 = ex1-1
        
        ## delete non-protein atoms
        self.ref = self.prepareRef( ref )

        ## delete non-protein atoms, cast to ref, fit to average and ref
        self.traj = self.prepareTraj( traj, self.ref, cast=self.cast )
        self.nframes = len( self.traj )

        ## reset frame limits for AmberCrdEntropist
        self.start = 0
        self.stop = len( self.traj )
        self.step = 1


    def prepare( self ):
        """Overrides Executor method."""

        ## create Parm file for selected atoms of trajectory
        if not os.path.exists( self.f_parm ):
            self.buildParm()
        else:
            self.log.add('using existing %s' % self.f_parm)

        ## create Amber Crd file for ptraj
        if not os.path.exists( self.f_crd ):
            self.log.add_nobreak('Writing amber crd file...')
            self.traj.writeCrd( self.f_crd )
            self.log.add('done')
            ## release memory
            if self.memsave: self.traj = None
        else:
            self.log.add('using existing %s' % self.f_crd)
        

    def buildParm( self ):
        a = AmberParmBuilder( self.traj.ref, verbose=self.verbose,
                              debug=self.debug )
        self.log.add_nobreak('Building amber topology...')

        a.parmMirror( self.f_parm, self.parmcrd )
        self.log.add('Topology built')


    def fit( self, traj, refModel=None, mask=None, conv=1e-6 ):
        """
        Fit trajectory until convergence onto it's own average and then
        transform the average of all frames onto the reference.
        """
        self.fit_e = self.fit_e or len( traj )
        self.fit_s = self.fit_s or 0
        
        traj.fit( ref=traj.ref, mask=mask )

        m_avg = traj[self.fit_s : self.fit_e ].avgModel()

        ## fit on average until it's not getting better
        d = 1.
        dd= 1.
        while dd >= conv:

            traj.fit( ref=m_avg, mask=mask )
            m_new_avg = traj[self.fit_s : self.fit_e].avgModel()

            oldD, d    = d, m_avg.rms( m_new_avg, mask=mask )
            self.log.add( "rms difference: %f" % d )

            dd = oldD - d
            m_avg = m_new_avg
            
        ## transform trajectory en block onto reference
        if refModel:
            self.log.add('fitting trajectory en-block onto reference...')

            if refModel.atomNames() != traj.ref.atomNames():
                self.log.add('casting ref for fitting...')
                ref_i, i = refModel.compareAtoms( m_avg )

                refModel = refModel.take( ref_i )
                m_avg    = m_avg.take( i )

                if mask:   mask = N.take( mask, i )
                
            r, t = m_avg.transformation( refModel, mask )
            traj.transform( r, t )
            

    def __splitFilenames( self, f ):
        """split(traj1.dat+traj2.dat) -> (traj1.dat, traj2.dat)
        """
        if f.find("+") != -1 :
            split = f.find("+")
            f1 = f[:split]
            f2 = f[split+1:]

            return f1, f2
        return None

    def __cleanAtoms( self, m ):
        """Remove non protein atoms and H if needed."""
        m.keep( N.nonzero( m.maskProtein() ) )
        if self.heavy:
            m.keep( N.nonzero( m.maskHeavy() ) )
        return m


    def __getModel( self, f ):
        """
        Load PDBModel directly or extract it from Trajectory.
        f - str, file name of PDB file, pickled PDBModel, or Trajectory
        -> PDBModel
        !! IOError, if file does not exist
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

        raise EntropistError, 'unknown reference type'


    def prepareRef( self, fname ):

        if not fname:
            return None

        if self.__splitFilenames( fname ):
            f1, f2 = self.__splitFilenames( fname )
            m1, m2 = PDBModel( self.__getModel(f1) ), \
                     PDBModel( self.__getModel(f2) )

            ref = Complex( m1, m2 )
        else:
            ref = t.Load( fname )

        if isinstance( ref, Trajectory ):
            ref = ref.ref

        if isinstance( ref, PDBModel ):
            return self.__cleanAtoms( ref )

        if isinstance( ref, Complex ):
            self.__cleanAtoms( ref.rec_model )
            self.__cleanAtoms( ref.lig_model )
            ref.lig_model_transformed = None
            return ref

        raise EntropistError, 'unknown reference type'


    def __add3( self, n_members, excluded, trippleIndex ):
        """
        Add a tripple of numbers from range( n_members ) to excluded.
        """
        remaining = MU.difference( range( n_members ), excluded )
        tripple = self.tripples( remaining, trippleIndex+1 )[-1]
        return MU.union( excluded, list(tripple) )

    def __add1( self, n_members, excluded, index ):
        """
        Add one number from range( n_members ) to list of excluded indices
        """
        remaining = MU.difference( range( n_members ), excluded )
        new_i = remaining[index]
        return excluded + [ new_i ]

    def __exclude( self, traj, exclude ):
        if exclude == None or len( exclude ) == 0:
            return traj

        members = range( traj.n_members )
        self.log.add("excluding members: " + str(exclude))

        return traj.takeMembers( MU.difference( members, exclude ) )
            
    def __removeMembers( self, t ):
        """
        t - EnsembleTraj OR ( EnsembleTraj, EnsembleTraj )
        """
        if self.ex_n:
            self.exclude = range( self.ex_n )

        if type( t ) is tuple and not type( self.exclude ) is tuple:
            self.exclude = ( self.exclude, self.exclude )

        if type( t ) is tuple:
            n_memb = ( t[0].n_members, t[1].n_members )
        else:
            n_memb = t.n_members
            
        if self.ex3 != None and type( self.exclude ) is tuple:
            self.exclude = self.__add3(n_memb[0], self.exclude[0], self.ex3),\
                           self.__add3(n_memb[1], self.exclude[1], self.ex3)

        if self.ex3 != None and type( self.exclude ) is list:
            self.exclude = self.__add3( n_memb, self.exclude, self.ex3 )

        if self.ex1 != None and type( self.exclude ) is tuple:
            self.exclude = self.__add1(n_memb[0], self.exclude[0], self.ex1),\
                           self.__add1(n_memb[1], self.exclude[1], self.ex1)

        if self.ex1 != None and type( self.exclude ) is list:
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
            raise EntropistError,'split,shift,shuffle require -border.'

        ## adapt reference to type of trajectory input 
        if ref and type(t) is tuple and not isinstance( ref, Complex ):
            rec = ref.takeChains( range(t[0].lenChains()) )
            lig = ref.takeChains( range(t[0].lenChains(), t[1].lenChains()) )
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

        self.log.add( 'Analysing trajectory with %i atoms and %i frames.' \
                      % (t.lenAtoms(), t.lenFrames()))

        return t


    def load_locked( self, fname ):
        """
        wait with unpickling until another Entropist has finished.
        """
        flock = fname + '__locked'

        while os.path.exists( flock ):
            self.log.add_nobreak('~')
            time.sleep( random.random() * 10 )
        self.log.add('')

        try:
            f = open( flock, 'w' )
            f.write('1')
            f.close()

            r = t.Load(fname)

        finally:
            t.tryRemove( flock )

        return r


    def loadTraj( self, fname, shift=0 ):
        """Load single trajectory."""
        if self.verbose:
            self.log.add_nobreak( 'Loading %s...' % fname)
            
        traj = self.load_locked( fname )

        if self.verbose:
            self.log.add( 'Processing trajectory...')

        ## convert single member frame index into all member frame index
        if (self.sstart or self.sstop) and isinstance(traj, EnsembleTraj):
            self.start = (self.sstart or 0) * traj.n_members
            self.stop  = (self.sstop  or 0) * traj.n_members
        if (self.sstart or self.sstop) and not isinstance(traj, EnsembleTraj):
            self.start, self.stop = self.sstart, self.sstop
            self.log.add('Warning: I am using -ss -se instead of -s -e')
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
            self.log.add( "Thinned to %i frames." % len( traj ) )
            
        ## keep only allowed atoms (default: all)
        if self.atoms:
            traj.ref.addChainId()
            aMask = traj.ref.mask( lambda a,ok=self.atoms: a['name'] in ok )
            traj.removeAtoms( N.nonzero( N.logical_not( aMask ) )  )

        ## get rid of non-standard atoms, water, ions, etc.
        l = traj.lenAtoms()
        traj = traj.compressAtoms( traj.ref.maskProtein() )
        self.log.add('%i non-protein atoms deleted.' % (l - traj.lenAtoms()) )

        ## delete hydrogens, if requested
        if self.heavy:
            l = traj.lenAtoms()
            traj = traj.compressAtoms( traj.ref.maskHeavy() )
            self.log.add('%i hydrogens deleted.' % (l - traj.lenAtoms()) )

        return traj


    def shiftTraj( self, traj, shift=0 ):
        """reorder member trajectories"""
        if not shift:
            return traj

        if not isinstance( traj, EnsembleTraj):
            raise EntropistError, 'shift requires EnsembleTraj'

        r = range( shift, traj.n_members ) + range( shift )
        self.log.add('reorder member trajectories: %s...' % str( r ) )

        return traj.takeMembers( r )


    def shuffleTraj( self, traj ):
        """reorder all frames at random"""
        r = range( len( traj ) )
        random.shuffle( r )
        return traj.takeFrames( r )

    def castTraj( self, traj, refModel ):
        """Equalize atom content of traj to refModel. """
        l    = traj.lenAtoms()
        lref = len( refModel )

        self.log.add('comparing traj with %i atoms to reference of %i atoms.'%\
                     (l, lref))
        i, iRef = traj.ref.compareAtoms( refModel )

        refModel.keep( iRef )
        self.log.add("%i atoms deleted from reference."%(lref-len(refModel)))

        traj.keepAtoms( i )
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
        pw = N.zeros( (len(all), len(all)), 'f' )
        for i in range( len( all ) ):
            for j in range( i, len(all) ):
                pw[i,j] = pw[j,i] = len( MU.intersection(all[i],all[j]) )**2

        pos = 0
        r = []

        while len( r ) < n:

            r += [ pos ]
            ## overlap of selected tripples with all others
            overlap = N.sum( N.array( [ pw[ i ] for i in r ] ) )
            ## select one with lowest overlap to all tripples selected before
            pos = N.argmin( overlap )

        return N.take( all, r )


    def cleanup( self ):
        AmberCrdEntropist.cleanup( self )

        if not self.debug:
            t.tryRemove( self.parmcrd )

            if not self.keep_crd:
                t.tryRemove( self.f_crd )
            if not self.keep_parm:
                t.tryRemove( self.f_parm )


    def finish( self ):
        AmberCrdEntropist.finish( self )

        ## consistency check of result
        if self.result['nframes'] != self.nframes:
            raise EntropistError, 'incorrect number of frames: %i instead %i'%\
                  ( self.result['nframes'], self.nframes )