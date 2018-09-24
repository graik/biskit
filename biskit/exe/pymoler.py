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

"""
Display Structures with Pymol
"""
## allow relative imports when calling module by itself for testing (pep-0366)
if __name__ == "__main__" and __package__ is None:
    import biskit.exe; __package__ = "biskit.exe"

import os
import tempfile

import biskit.tools as T
import biskit.settings as settings
from biskit import PDBModel

from .executor import Executor, TemplateError


## ================= single Pymol model ===========================

class PymolModel:

    def __init__( self, model, modName ):
        """
        @param model: model to view
        @type  model: PDBModel
        @param modName: model name, will show up in PyMol
        @type  modName: str
        """
        self.fname = ''
        self.temporary = 0
        self.struct = None

        if type( model ) is str:
            self.fname = model
        else:
            self.struct = model.clone() ## clone atom dicts
            self.temporary = 1

        self.modName = modName

        if self.fname == '':
            self.fname = tempfile.mktemp( 'model')


    def addProperty( self, values, key='temperature_factor' ):
        """
        Add extra value to each atom in Structure. The values
        will be written to either the B- (temperature_factor)
        or Q-factor 'occupancy' column in the temporary pdb-file.
        These values can then be used to display properties in PyMol
        via commands like 'color_b' and 'color_q'. See also
        L{addResProperty}.

        @param values: list of numbers, len( values ) == number of atoms
        @type  values: [float]      
        @param key: key for Atom.properties dictionary
                    ('occupancy' OR 'temperature_factor')
        @type  key: occupancy|temperature_factor
        """
        if self.struct is None:
            self.struct = PDBModel( self.fname )
            self.temporary = 1

        self.struct[ key ] = values


    def addResProperty( self, values, key='temperature_factor'):
        """
        Does the same thing as L{addProperty} but on the residue level,
        i.e adds extra value to each residue in Structure.
        (The same value is added to all atoms of a residue.)
        These values can then be used to display properties in PyMol
        via commands like 'color_b' and 'color_q'.
        
        @param values: list of numbers, len( values ) == number of residues
        @type  values: [float]      
        @param key: key for Atom.properties dictionary
                    ('occupancy' OR 'temperature_factor')
        @type  key: occupancy|temperature_factor
        """
        try:
            if self.struct is None:
                self.struct = PDBModel( self.fname )
                self.temporary = 1

            self.struct[ key ] = self.struct.res2atomProfile( values )
        except:
            print(T.lastError())


    def writeIfNeeded( self ):
        """
        Create pdb on disc if it's not already there or if it
        has to be changed.
        
        @return: filename of new or existing pdb
        @rtype: str
        """
        if self.temporary:
            self.struct.writePdb( self.fname, 2 )
        return self.fname


## ========================= Pymoler ============================
## create Pymol input script and display in Pymol ##

class Pymoler( Executor ):
    """
    Run Pymol
    =========
    Create a pymol script file (.pml) and execule it using Pymol.

    Example usage
    -------------
        >>> pm = Pymoler()
        >>> pm.addPdb( model )
        >>> pm.run()

    References
    ----------
       - U{http://pymol.sourceforge.net/}
    
    @note: the file is only flushed to disc when the object
           is destructed, or the flush() method is called!
    @note: since version 2.6 Pymoler is using Biskit.Executor
    """

    def __init__(self, full=0, mode='w', verbose=1, **kw ):
        """
        @param mode: open file with this mode, w=override, a=append
        @type  mode: str
        @param full: dispaly pymol structures in fill screen mode::
                       0 - normal mode
                       1 - full screen mode
                       2 - full screen and no menues
        @type  full: 0|1|2        
        """
        self.verbose = verbose
        
        # name of .pml file
        self.foutName = tempfile.mktemp() + '.pml'

        # open for <appending|writing|reading>
        self.fgenerate = open(self.foutName, mode) 

        # will contain PymolModels or lists of PymolModels
        self.dic = {}

        ## add startup commands
        self.initPymol()

        ## set arguments for display options (normal, full, all)
        arg = '-q %s'%self.foutName
        if full == 1:
            arg = '-qe %s'%self.foutName
        if full == 2:
            arg = '-qei %s'%self.foutName
            
        Executor.__init__( self, 'pymol', args=arg,
                           catch_err=1, catch_out=1, **kw )
        
    def __del__(self):
        self.fgenerate.close()


    def flush(self):
        """
        Flush output file (but keep it open).
        """
        self.fgenerate.flush()


    def add(self, str):
        """
        Add String str and line break to file.

        @param str: string to add to pml file
        @type  str: str        
        """
        try:
            self.fgenerate.write(str + '\n')
        except (IOError):
            T.errWriteln(
                "PymolInput.add(): Error adding string to pymol script file.")
            T.errWriteln( T.lastError() )


    def addMovie( self, pdb, modName=None ):
        """
        Add one or several existing pdb files or Structure objects
        to one model. Several files will hence end up as single movie
        (i.e. as frames of a model in PyMol).
 
        @param pdb: file name or a list of file names OR
                    PDBModel or list of PDBModels
        @type  pdb: str or [str] OR PDBModel or [PDBModel]
        @param modName: model name, will show up in PyMol. If 'None' a
                        model name will be created from the source file
                        name and a serial number.
        @type  modName: str OR None

        @return: the modName of the added model
        @rtype: str        
        """
        if type( pdb ) is not list:
            pdb = [pdb]

        ## dream up a nice model name
        if modName is None:

            if type( pdb[0]) is str:
                modName = T.stripFilename( pdb[0] )
                modName = self._getFreeModName( modName, 0 )

            else:
                modName = self._getFreeModName( 'models', 0 )

        ## create new empty list
        if modName not in self.dic:
            self.dic[ modName ] = []

        ## create model object for each file and append it to dic
        for f in pdb:

            ## create model from Structure or file name
            model = PymolModel( f, modName )
            self.dic[ modName ].append( model )

            ## add load statement to Pymol script
            self.add( 'load '+ model.fname + ',' + modName )

        return modName


    def addFrame( self, frame, modName ):
        """
        Add file(s) or Structure(s) to an EXISTING model as movie.
        
        @param frame: the structure to add to an existing model
        @type  frame: str(s) or PDBModel(s) 
        @param modName: model name, must be existing
                        (i.e. an already added model)
        @type  modName: str
        """
        self.addMovie( frame, modName )


    def _getFreeModName( self, base, index ):
        """
        Return next free model name in dictionary, made up
        from base + index
        
        @param base: name base
        @type  base: str
        @param index: name suffix 
        @type  index: int
        
        @return: name composed of base+suffix
        @rtype: str
        """
        if base in self.dic:

            if base + str(index) in self.dic:

                return self._getFreeModName( base, index + 1 )

            else:
                return base + str( index )

        else:
            return base


    def addPdb( self, pdb, modName=None ):
        """
        Add one or several existing pdbs. Make sure all go into
        different models (have different model names).
        
        @param pdb: file name or a list of file names OR
                    PDBModel or list of PDBModels
        @type  pdb: str or [str] OR PDBModel or [PDBModel]

        @param modName: force model name, will change for list
                        of file names. If 'None' a model name will
                        be created from the source file name.
        @type  modName: str OR None
        
        @return: model name of first file
        @rtype: str
        """

        if type( pdb ) is not list:
            pdb = [pdb]

        ## create one model for each file / Structure
        for f in pdb:
            result = self.addMovie( f, modName )

        return result


    def makeSel( self, selDic ):
        """
        Make a selection.
        
        @param selDic: a selection dictionary,  that can be of three types:
                         1. dictionary with one or more of the keys:
                           'model' 'segment', 'chain', 'residue', 'atom'
                         2. dictionary with key:   'element'
                         3. dictionaty with key:   'expression'
        @type  selDic: dict
        """
        ## the selection must comply with one of the three selection
        ## schemes below
        keyList = ['model', 'segment', 'chain', 'residue', 'atom', \
                   'element', 'expression']

        # check for invalid selections
        for sel in selDic.keys():
            if sel not in keyList:
                print('Invalid selection in ' + str(selDic))

        # element
        if keyList[5] in selDic:
            selection = '( elem ' + str(selDic[keyList[5]]) + ' )'

        # expression
        elif keyList[6] in selDic:
            selection = '( ' + str(selDic[keyList[6]]) + ' )'

        # create selection from selDic
        else:
            for key in keyList:
                if key not in selDic:
                    selDic[key]=''
            selection = '( /' + str(selDic['model']) + \
                        '/' + str(selDic['segment']) + \
                        '/' + str(selDic['chain']) + \
                        '/' + str(selDic['residue']) + \
                        '/' + str(selDic['atom']) + ' )'

        return selection


    def writeStructures( self ):
        """
        Write all needed PDB files to disc.
        """
        for key in self.dic.keys():

            for model in self.dic[ key ]:

                n = model.writeIfNeeded()
                
                if self.verbose: print(n)


    def deletePdbs(self):
        """
        Deletes the pdb-files in the list from disc
        """
        for key in self.dic.keys():

            for model in self.dic[ key ]:
                ## only remove files created by PymolInput!
                if model.temporary:
                    T.tryRemove(model.fname)


    def setAtomValues( self, model, values, key='temperature_factor',
                       lastOnly=0 ):
        """
        Add numeric value to all atoms of all Structures or the last
        Structure in a model.. The values will be written to either
        the B- (temperature_factor) or Q-factor 'occupancy' column in
        the temporary pdb-file.
        These values can then be used to display properties in PyMol
        via commands like 'color_b' and 'color_q'. See also
        L{setResValues}.

        @param model: model name
        @type  model: str
        @param values: list of numbers, len( values ) == number of atoms
        @type  values: [float]      
        @param key: key for Atom.properties dictionary
                    (default: temperature_factor)
        @type  key: occupancy|temperature_factor
        @param lastOnly: 0 - add to all in model OR
                         1 - add only to last Structure (default: 0)
        @type  lastOnly: 1|0
        """
        if lastOnly:
            self.dic[ model ][-1].addProperty( values, key )

        else:
            for m in self.dic[ model ]:
                try:
                    m.addProperty( values, key )
                except:
                    T.errWriteln( "Warning: error while adding properties.")
                    T.errWriteln( "Key: "+str( key )+" values: "+str( values ) )
                    T.errWriteln( T.lastError() )


    def setResValues( self, model, values, key='temperature_factor',
                      lastOnly=0 ):
        """
        Add numeric value per residue to all atoms of all Structures
        or the last Structure in a model. The values will be written to
        either the B- (temperature_factor) or Q-factor 'occupancy' column
        in the temporary pdb-file.
        These values can then be used to display properties in PyMol
        via commands like 'color_b' and 'color_q'. See also
        L{setAtomValues}.
        
        @param model: model name
        @type  model: str
        @param values: list of numbers, len( values ) == number of residues
        @type  values: [float]      
        @param key: key for Atom.properties dictionary
                    (default: temperature_factor)
        @type  key: occupancy|temperature_factor
        @param lastOnly: 0 - add to all in model OR
                         1 - add only to last Structure (default: 0)
        @type  lastOnly: 1|0
        """
        if lastOnly:
            self.dic[ model ][-1].addResProperty( values, key )

        else:
            for m in self.dic[ model ]:
                try:
                    m.addResProperty( values, key )
                except:
                    T.errWriteln( "Warning: error while adding properties.")
                    T.errWriteln( "Key: "+str( key )+" values: "+str( values ) )
                    T.errWriteln( T.lastError() )


    def colorAtoms( self, model, values, lastOnly=0 ):
        """
        Color atoms of this model by list of values.
        
        @param model: model name
        @type  model: str
        @param values: len == number of atoms
        @type  values: [float]
        @param lastOnly: 0 - add to all in model OR
                         1 - add only to last Structure (default: 0)
        @type  lastOnly: 1|0         
        """
        self.setAtomValues( model, values, lastOnly=lastOnly )
        self.add( "color_b('%s')" % model )


    def colorRes( self, model, values, lastOnly=0 ):
        """
        Color residues of this model by list of values.
        
        @param model: model name
        @type  model: str
        @param values: len == number of residues
        @type  values: list of numbers
        @param lastOnly: 0 .. add to all in model
                         1 .. add only to last Structure
        @type  lastOnly: 
        """
        self.setResValues( model, values, lastOnly=lastOnly )
        self.add( "color_b('%s')" % model )



    def addColors( self, nColors, firstColor=[1.0, 0.0, 0.0],
                   lastColor=[0.0, 1.0, 0.0] ):
        """
        Define a range of colors that can be called in PyMol by name.
        
        @param nColors: numbers of colors generated
        @type  nColors: int
        @param firstColor: first rgb color (default: [1.0, 0.0, 0.0])
        @type  firstColor: [float]
        @param lastColor: last rgb color (default: [0.0, 1.0, 0.0])
        @type  lastColor: [float]
        @return: a list of the color names
        @rtype: [str]
        """   
        spectrum = T.hexColors( nColors, T.rgb2hex( firstColor ),
                              T.rgb2hex( lastColor ) )
        rgb = []
        colorNames = []

        for c in range(0, nColors):
            rgb = [ float( T.hex2int(spectrum[c][:4]) ) / 255, \
                    float( T.hex2int('0x' + spectrum[c][4:6]) ) / 255, \
                    float( T.hex2int('0x' + spectrum[c][6:]) ) / 255 ]

            cName = 'c' + str(c)
            colorNames += cName

            self.add( 'set_color '+ cName + ', ' + str(rgb) )

        return colorNames


    def initPymol( self ):
        """
        Do some stuff always first...
        """
        ## import dssp command and other commands
        if os.path.isdir( settings.pymol_scripts ):
            for script in os.listdir( settings.pymol_scripts ):
                if not os.path.isdir( script ):
                    self.add( 'run ' + settings.pymol_scripts + script )
                    if self.verbose:
                        print('Adding %s script as a PyMol command'%script)
        else:
            print('\n\nWARNING: No external PyMol scripts added\n\n')


##      # MAKE default SELECTIONS
##      self.add('select bb, ' + self.makeSel({'expression':'*/c,ca,o,n'} ))
##      self.add('select sc, ' +\
##               self.makeSel({'expression':'!(*/c,o,ca,n) and !element h'}))
##         self.add('select none')


    def prepare( self, cleanUp=1 ):
        """
        Overrides Executor method.
        """
        Executor.prepare(self)
        
        self.flush() ## important to avoid empty input file

        ## Write PDB's to disc if needed
        self.writeStructures()

    def cleanup(self):
        if not self.debug:
            T.tryRemove( self.foutName )
            self.deletePdbs()
    
        Executor.cleanup( self )
        

    def show( self ):
        """
        Backward compatability with old scripts.
        """
        self.run()


#############
##  TESTING        
#############
import biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test"""

    TAGS = [ BT.EXE ]

    def test_Pymoler(self):
        """Pymoler test"""
        self.traj = T.load( T.testRoot() + '/lig_pcr_00/traj.dat' )

        self.pm = Pymoler( full=0, verbose=self.local )
        
        mname = self.pm.addMovie( [ self.traj[i] for i in range(0,100,20) ] )

        sel = self.pm.makeSel({'residue':29})
##         self.pm.add('show stick, %s' % sel)
        self.pm.add('show surface, %s' % sel)

        self.pm.add('mplay')
        
        if not self.local:
            self.pm.add('quit')
            
        self.pm.run() ## old style call "pm.show()" also works
    
        self.assertTrue( self.pm.pid is not None )
        
if __name__ == '__main__':

    BT.localTest()



