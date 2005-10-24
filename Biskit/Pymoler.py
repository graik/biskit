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
## last $Date$
## $Revision$
"""Display Structures with Pymol"""


from PDBModel import PDBModel
import tools as T
import os
import tempfile
import settings

# =============== SETUP PYMOL ENVIRONMENT =========================

## check the pymol.com file for path settings

## # dynamic linking 
## try:
##     os.environ['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH'] + ':' \
## 					+ os.environ['PYMOL_EXTLIBPATH']
## except:
##     os.environ['LD_LIBRARY_PATH'] = os.environ['PYMOL_EXTLIBPATH']
	
## # python modules
## try: 
##     os.environ['PYTHONPATH'] = os.environ['PYMOL_PATH'] + ':' \
## 				   + os.environ['PYMOL_EXTLIBPATH'] \
## 				   + 'python2.2' \
## 				   + ':' +os.environ['PYTHONPATH']
## except:
##     os.environ['PYTHONPATH'] = os.environ['PYMOL_PATH'] + ':' \
## 				   + os.environ['PYMOL_EXTLIBPATH'] \
## 				   + 'python2.2'        
## # temp directory path
os.environ['TMPDIR'] = settings.tempDirShared


## ================= single Pymol model ===========================

class PymolModel:

    def __init__( self, model, modName ):

        self.fname = ''
        self.temporary = 0
        self.struct = None
        
        if type( model ) is str:
            self.fname = model
        else:
            self.struct = model.clone( deepcopy=1 ) ## clone atom dicts
            self.temporary = 1

        self.modName = modName

        if self.fname == '':
            self.fname = tempfile.mktemp( 'model')


    def addProperty( self, values, key='temperature_factor' ):
        """
        Add extra value to each atom in Structure.
        key - string, key for Atom.properties dictionary
	      ('occupancy', 'temperature_factor')
        values - list of numbers, len( values ) == atom number
        """
        if self.struct == None:
            self.struct = PDBModel( self.fname )
            self.temporary = 1

        i = 0
        for atom in self.struct.atoms:
            atom[ key ] = values[i]
            i += 1

                
    def addResProperty( self, values, key='temperature_factor'):
        """
        Add extra value to each residue in Structure.
        (The same value is added to all atoms of a residue.)
        key - string, key for Atom.properties dictionary
	      ('occupancy', 'temperature_factor')
        values - list of numbers, len( values ) == residue number
        """
        try:
            if self.struct == None:
                self.struct = PDBModel( self.fname )
                self.temporary = 1

            i = 0
            resMap = self.struct.resMap()

            for atom in self.struct.atoms:
                atom[ key ] = values[ resMap[i] ]
                i += 1
        except:
            print T.lastError()


    def writeIfNeeded( self ):
        """
        Create pdb on disc if it's not already there or if it
        has to be changed.
        -> string, filename of new or existing pdb
        """
        if self.temporary:
            self.struct.writePdb( self.fname, 2 )
        return self.fname

## ========================= Pymoler ============================
## create Pymol input script and display Pymol ##

class Pymoler:
    """
    Create PyMol script (.pml) file.
    Note, the file is only flushed to disc when the object is destructed,
    or the flush() method is called!
    """

    def __init__(self, mode='w'):
        """
        mode    - open file with this mode, w=override, a=append
        """
        #self.foutName = os.path.abspath(outName)    # name of new .pml file
	self.foutName = tempfile.mktemp() + '.pml'
	
        # open for <appending|writing|reading>
        self.fgenerate = open(self.foutName, mode) 

        # will contain PymolModels or lists of PymolModels
        self.dic = {}

        ## add startup commands
        self.initPymol()


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
        to one model. Several files will hence end up as single movie.
        
        pdb - str or list of str, file name(s) or ..
            - PDBModel or list of it
        modName - str, model name
        """
        if type( pdb ) is not list:
            pdb = [pdb]

        ## dream up a nice model name
        if modName == None:

            if type( pdb[0]) is str:
                modName = T.stripFilename( pdb[0] )
                modName = self._getFreeModName( modName, 0 )

            else:
                modName = self._getFreeModName( 'model', 0 )

        ## create new empty list
        if not self.dic.has_key( modName ):
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
        Add file(s) or Structure(s) to EXISTING model as movie.
        frame - str or Scientific.IO.PDB.Structure
        modName - str, model name, must be existing
        """
        self.addMovie( frame, modName )


    def _getFreeModName( self, base, index ):
        """
        Return next free model name in dictionary, made up
        from base + index
        base - str
        index - int
        -> str
        """
        if self.dic.has_key( base ):

            if self.dic.has_key( base + str(index) ):

                return self._getFreeModName( base, index + 1 )

            else:
                return base + str( index )

        else:
            return base


    def addPdb( self, pdb, modName=None ):
        """
        Add one or several existing pdbs. Make sure all go into
        different model (have different model name).
        pdb - str or [str, str, ..], file name or list of names
              - PDBModel or list of it
        modName - str, force model name, will change for list
                  of file names
        -> str, model name of first file
        """

        if type( pdb ) is not list:
            pdb = [pdb]

        ## create one model for each file / Structure
        for f in pdb:
            result = self.addMovie( f, modName )

        return result
    

    def makeSel( self, selDic ):
        """
        Make a selection, that can be of three types:
           1) dictionary with one or more of the keys:
                'model' 'segment', 'chain', 'residue', 'atom'
           2) dictionary with key:   'element'
           3) dictionaty with key:   'expression'
        """
        ## the selection must comply with one of the three selection
	## schemes below
        keyList = ['model', 'segment', 'chain', 'residue', 'atom', \
		   'element', 'expression']

        # check for invalid selections
        for sel in selDic.keys():
            if sel not in keyList:
                print 'Invalid selection in ' + str(selDic)
        
        # element
        if selDic.has_key(keyList[5]):
            selection = '( elem ' + str(selDic[keyList[5]]) + ' )'

        # expression
        elif selDic.has_key(keyList[6]):
            selection = '( ' + str(selDic[keyList[6]]) + ' )'

        # create selection from selDic
        else:
            for key in keyList:
                if not selDic.has_key(key):
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
                print n

    def addDeleteScript(self):
        """
	Deletes the pymol script file from disc
	"""
	self.add( "/ import os" )
	self.add("/ os.system('rm " + self.foutName+ "')")


    def addDeletePdbs(self):
        """
	Deletes the pdb-files in the list from disc
	"""
	self.add( "/ import os" )
	for key in self.dic.keys():

            for model in self.dic[ key ]:
                ## only remove files created by PymolInput!
		if model.temporary:
		        self.add( "/ os.system('rm " + model.fname + "')" )


    def setAtomValues( self, model, values, key='temperature_factor',
		       lastOnly=0 ):
        """
        Add numeric value to all atoms of all Structures
        or the last Structure in a model.
        model - str, model name
        values - list of numbers, len == number of atoms
        key - key for atom.properties dictionary
        lastOnly - 0 .. add to all in model
                   1 .. add only to last Structure
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
        or the last Structure in a model.
        model - str, model name
        values - list of numbers, len == number of residues
        key - key for atom.properties dictionary
        lastOnly - 0 .. add to all in model
                   1 .. add only to last Structure
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
        model - str, model name
        values - list of numbers, len == number of atoms
        lastOnly - 0 .. add to all in model
                   1 .. add only to last Structure
        """
        self.setAtomValues( model, values, lastOnly=lastOnly )
        self.add( "color_b('%s')" % model )
        

    def colorRes( self, model, values, lastOnly=0 ):
        """
        Color residues of this model by list of values.
        model - str, model name
        values - list of numbers, len == number of residues
        lastOnly - 0 .. add to all in model
                   1 .. add only to last Structure
        """
        self.setResValues( model, values, lastOnly=lastOnly )
        self.add( "color_b('%s')" % model )
        


    def addColors( self, nColors, firstColor=[1.0, 0.0, 0.0],
                   lastColor=[0.0, 1.0, 0.0] ):
        """
        Define a range of colors to be called in PyMol by name.
        nColors - numbers of colors generated
        firstColor - first rgb color
        lastColor - last   ~    ~
        -> colorNames, a list of the color names
        """   
        spectrum = T.hexColors( nColors, T.rgb2hex( firstColor ),
                              T.rgb2hex( lastColor ) )
        rgb = []
        colorNames = []

        for c in range(0, nColors):
            rgb = [ float( spectrum[c][:4] ) / 255, \
                  float( '0x' + spectrum[c][4:6] ) / 255, \
                  float( '0x' + spectrum[c][6:] ) / 255 ]

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
                    print 'Adding %s script as a PyMol command'%script
                    self.add( 'run ' + settings.pymol_scripts + script )
        else:
            print '\n\nWARNING: No external PyMol scripts added\n\n'
    

##  	# MAKE default SELECTIONS
## 	self.add('select bb, ' + self.makeSel({'expression':'*/c,ca,o,n'} ))
##  	self.add('select sc, ' +\
## 		 self.makeSel({'expression':'!(*/c,o,ca,n) and !element h'}))
##         self.add('select none')
        
    
    def show(self, cleanUp=1 ):
	"""
	Takes a structure object.
	Dispays it in PyMol, files are only temporarily created.
	"""
	# clean up files from disc
        if cleanUp:
            self.addDeletePdbs()
            self.addDeleteScript()
	
	self.flush()

        ## Write PDB's to disc if needed
	self.writeStructures()
        
        ## -q for quiet mode
	os.spawnlpe(os.P_NOWAIT, settings.pymol_bin,
                    settings.pymol_bin , '-q', self.foutName, os.environ )


    def fshow(self, cleanUp=1, sideMenu=1 ):
	"""
	Dispays all added Pdbs full-screen in PyMol,
        files are only temporarily created.
        cleanUp - 0||1, remove temporary files
        sideMenu- 0||1, show side menu
	"""
	# clean up files from disc
        if cleanUp:
            self.addDeletePdbs()
            self.addDeleteScript()
	
	self.flush()

        ## Write PDB's to disc if needed
	self.writeStructures()

        options = '-qex'
        if not sideMenu:
            options += 'i'

  	os.spawnlpe(os.P_NOWAIT, settings.pymol_bin, settings.pymol_bin,
                    options, self.foutName, os.environ )

## TEST ##

if __name__ == '__main__':

    import glob
    
    f = glob.glob( T.testRoot()+'/lig_pc2_00/pdb/*_1_*pdb.gz' )[:2]
    m = [ PDBModel(i) for i in f ]
    m = [ i.compress( i.maskProtein() ) for i in m ]
    
    pm = Pymoler( )
    pm.addMovie( [ m[0], m[1] ] )

    pm.show()
