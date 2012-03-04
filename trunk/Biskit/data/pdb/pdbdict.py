import numpy as N
import Biskit.tools as T

class Residue( object ):
    """
    Description of the standard content of a certain residue.
    """

    #: residue type codes
    UNKNOWN = 0
    PROTEIN = 1
    DNA     = 2
    RNA     = 3
    HETATM  = 10

    def __init__(self, name, atoms,
                 single=None, type=None, elements=[], masses=[] ):
        """
        """
        self.name = name       #: typical residue code
        self.single = single   #: single-letter residue code
        self.type = type or self.UNKNOWN
        self.atoms = N.array( atoms, dtype=str )
        self.masses = N.array( masses, dtype=float ) 
        self.elements = N.array( elements, dtype=str ) 

    def _getmass( self ):
        return N.sum( self.masses )

    mass = property( _getmass )


class ResidueError(Exception):
    pass

class ResidueTranslator( object ):
    """
    Translate between PDB version 2.3 and 3.0 residue/atom names
    """

    FTABLE = T.dataRoot() + '/pdb/master_hash.txt'

    def __init__(self, ftable=FTABLE ):
        self.ftable = ftable
        self.resdic = {}  #: new to old translation of residue names
        self.atmdic = {}  #: new to old translation of atom names

        self.__parseTable()


    def __parseTable( self ):
        try:
            for line in open( self.ftable ).readlines():
                rnew, anew, rold, aold = self.__parseLine( line )

                if not rnew in self.resdic:
                    self.resdic[rnew] = rold

                if not rnew in self.atmdic:
                    self.atmdic[ rnew ] = {}

                self.atmdic[rnew][anew] = aold

        except IOError, why:
            raise ResidueError, 'Cannot read residue translation table %s'\
                  % self.ftable

            
    def __parseLine( self, l ):
        l = l.strip()
        new, old = l.split(':')
        aold, rold = old.split()
        anew, rnew = new.split()

        return rnew, anew, rold, aold

    def new2old( self, res ):
        """
        @type res: Residue
        @rtype: Residue
        """
        r = copy.copy( res )
        trans = self.atmdic[ r ]
        r.atoms = [ trans[a] for a in r.atoms ]
        r.name  = self.resdic[ r.name ]

        return r


rt = ResidueTranslator()

atoms = {}
for r, rdic in rt.atmdic.items():

    for a in rdic.keys():

        atoms[a] = rdic[a]
        
