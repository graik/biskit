from __future__ import with_statement
import os
import re
import sys
import tempfile
import subprocess

import Biskit.ProfileCollection as BP
import numpy as N

class LeapObjectFile:
    
    def __init__(self, file_name ):
        
        self.boxes = []
        self.molecules = []
        self.params = []
        self.all = {}
        
        try:
            file_H = open(file_name,'r')
            file_H.close()
        except IOError:
            sys.exit("File Error: %s"%(file_name))
            
        self.source = file_name
        self.parseAll()
    
    def __getitem__(self, name):
        """Direct access to unit instances in chosen list"""
        return self.all[name]
        
    def get(self, name):
        try:
            return self[name]
        except:
            return False
        
    def parseAll(self, file_name = None):
        
        if file_name is None:
            file_name = self.source
        
        UnitNames = self.getNames(file_name)
        [ self.getUnits(UnitName, file_name) for UnitName in UnitNames ]
        
    def getUnits(self, UnitName, file_name = None):
        
        if file_name is None:
            file_name = self.source
            
        OFF = open(file_name, 'r')
        search = "!entry."+UnitName+"."
        line = findNext( OFF, search, return_line = True )
        line = line.strip()
        OFF.close()
                
        Type = line.split('.')[2]
        if Type in 'unit':
            
            U_instance = Unit(UnitName, file_name)
            
            if U_instance.hasBox():
                U_instance.kind = 'Box'
                self.boxes.append( U_instance )
            
            else:
                U_instance.kind = 'Molecule'
                self.molecules.append( U_instance )
        
        elif Type in 'parm':
            U_instance = Parameter(UnitName, file_name)
            U_instance.kind = 'Parm'
            self.params.append( U_instance )
        
        self.all[U_instance.name]=U_instance
        
    def getNames(self, file_name = None):
        
        if file_name is None:
            file_name = self.source
        
        OFF = open(file_name, 'r')
        
        #Read Unit names
        OFF.readline()
        
        UnitNames = []
        while 1:
          line = OFF.readline().strip()
          if "!" in line: break
          U_name = (line.split("\""))[1]
          UnitNames.append(U_name)
        
        OFF.close()

        return UnitNames

class Unit:
    
    def __init__(self, name, file_name = None ):
        self.type = "Unit"
        self.name = name
        self.description = None
        self.residues = BP()
        self.res_num = None
        self.atoms = BP()
        self.box = {}
        
        if file_name is not None:
            #If file specified, read unit info
            self.read(file_name)
        
    def __repr__(self):
        return self.description
        
    def __getitem__( self, k ):
        """
        Get atom profile or profile item or CrossView for one atom::
          m['prof1']         <==>  m.atoms.get( 'prof1' )         
          m['prof1','info1'] <==>  m.atoms.get( 'prof1','info1' )
          m[10]              <==>  CrossView( m.atoms, 10 )
        
        @return: profile OR meta infos thereof OR CrossView dict
        @rtype: list OR array OR any OR CrossView
        """
        if type( k ) is str:
            if k in self.atoms:
                return self.atoms.get( k )
            if k in self.residues:
                return self.residues.get( k )
        
        return self.atoms[k]
    
    def hasBox(self):
        """Check if unit has box information"""
        return bool(self.box)
    
    def read(self, file_name):
        """
        Read all information for self in the Leap Object File file_name.
        Be aware, that we are not collecting atomspertinfo, childsequence, connect, connectivity, hierarchy, residueconnect, 
        solventcap nor velocities.
        We are also ignoring some fields in the explored blocks.
        Further developement is welcome if required.
        """
        OFF = open(file_name, 'r')
        self.source = file_name
        search = self.name + '.unit'
        
        ### ATOMS INFO
        # Head to atoms block
        findNext(OFF, search, return_line = False)
        
        atoms_block = []
        while 1:
            line = OFF.readline().strip()
            if search in line: break
            atoms_block.append(line.split(' '))
        atoms_block = N.array(atoms_block, dtype=N.str_)
        
        # Get rid of the double quotes in atom name and type
        for pos in (0,1):
          for i, atom in enumerate(atoms_block[:,pos]):
            atoms_block[i,pos] = atom.split('"')[1]
        
        self.atoms.set( 'name', list(atoms_block[:,0]) )
        self.atoms.set( 'atype', list(atoms_block[:,1]) )
        self.atoms.set( 'residue_number', list(atoms_block[:,3].astype(int)) )
        self.atoms.set( 'atom_id', list(atoms_block[:,5].astype(int)) )
        self.atoms.set( 'element', list(atoms_block[:,6].astype(int)) )
        self.atoms.set( 'charge', list(atoms_block[:,7].astype(float)) )
        
        
        
        ### BOX INFO
        # Skip atomspertinfo block and head to boundboux block
        findNext(OFF, search, return_line = False)
        
        _box = float(OFF.readline().strip())
        
        if _box == 1:
            angle = float(OFF.readline().strip())
            x = float(OFF.readline().strip())
            y = float(OFF.readline().strip())
            z = float(OFF.readline().strip())
            
            self.box['size'] = [x, y, z]
            self.box['volume'] = x*y*z
            self.box['angle'] = angle
        
        else:
            self.box = None
        
        ### DESCRIPTION
        # Skip 4 blocks and head to name block
        [ findNext(OFF, search, return_line = False) for it in xrange(5)]
        self.description = (OFF.readline().strip()).split('"')[1]
        
        
        ### COORDINATES
        findNext(OFF, search, return_line = False)
        
        positions_block = []
        while 1:
            line = OFF.readline().strip()
            if search in line: break
            positions_block.append(line.split(' '))
        positions_block = N.array(positions_block, dtype=float)
        
        self.atoms.set( 'coordinates', list(positions_block) )
        
        # Center of mass
        self.mcent = ( self.atoms['coordinates'].max(axis=0) - self.atoms['coordinates'].min(axis=0) ) / 2.
        
        
        ### RESIDUES
        #Skip next block and head to residues one
        findNext(OFF, search, return_line = False)
        
        residues_block = []
        while 1:
            line = OFF.readline().strip()
            if search in line: break
            residues_block.append(line.split(' '))
        residues_block = N.array(residues_block, dtype=str)
        
        #Get rid of the double quotes in residue name and type
        for pos in (0, 4):
          for i, res in enumerate(residues_block[:,pos]):
            residues_block[i,pos] = res.split('"')[1]
        
        self.residues.set( 'name', list( residues_block[:,0] ) )
        self.residues.set( 'id', list( residues_block[:,1].astype(int) ) )
        self.residues.set( 'head_atom', list( residues_block[:,3].astype(int) ) )
        self.residues.set( 'type', list( residues_block[:,4] ) )
        
        #Read Pdb Sequence numeration
        Pdb_resID = []
        while 1:
            line = OFF.readline().strip()
            if search in line: break
            Pdb_resID.append( int(line) )
        
        self.residues.set( 'PdbID', Pdb_resID )
        self.res_num = len(Pdb_resID)
        
        OFF.close()


class Parameter:
    
    def __init__(self, name, file_name = None):
        self.type = "Parameter"
        self.name = name
        self.description = name
        self.params = {}
        
        if file_name is not None:
            self.read(file_name)
    
    def __repr__(self):
        return self.description
    
    def read(self, file_name):
        """
        Read parameters Object.
        These blocks are more flexible than normal unit's.
        Can contain or not angle, bonds, atoms and/or torsions parameters.
        Information is stored but not sorted or treated.
        THIS NEEDS FURTHER DEVELOPEMENT.
        """
        search = self.name + '.parm'
        
        
        #First give a first round to identify number of entries
        entries = 0
        with open(file_name, 'r') as f:
            for line in f:
                if search in line:
                    entries += 1
        
        OFF = open(file_name, 'r')
        
        line = findNext(OFF, search, return_line = True)
        for i in range(entries):
            
            parm = line.split(' ')[0].split('.')[-1]
            block_info = []
            
            while 1:
              line = OFF.readline().strip()
              if re.match('!',line) or line == '': break
              block_info.append(line.split())
            
            self.params[parm] = block_info
                          
        OFF.close()


def findNext(file_handle, expr, return_line=True):
    
    match = None
    while 1:
        line = file_handle.readline()
        if expr in line: 
          match = 1
          break
        
    if match is None:
        file_handle.close()
        sys.exit("Couldn't find Unit information in the given Object File!")
    
    if return_line is False:
        return 1
    else:
        return line.strip()


class Leaper:
    """Handle Leap I/O"""
    def __init__(self,  objectFileName=None,  cwd=None):
        self.tleap = subprocess.Popen("tleap -f - ",  shell=True,   stdin = subprocess.PIPE,  stdout=subprocess.PIPE,  stderr=subprocess.PIPE, cwd=cwd)
        self._in = self.tleap.stdin
        self._out = self.tleap.stdout
        self._err = self.tleap.stderr
        if objectFileName: self.command('loadOff %s'%(objectFileName))

 
    def command(self,  command):
        """
        Command must be a string corresponding to a correct Leap command.
        Returns the output of leap for the given command.
        tLeap is a bit tricky to be controlled by stdin, it terminates when it reaches EOF.
        Thus we use a fake command to detect the last line before EOF.
        """
        self._in.write(command+'\n')
        self._in.write('lastCommandOut\n')
        out = []
        while 1:
            line = self._out.readline().strip()
            if 'ERROR: syntax error' in line: break
            if line: out.append(line)
        return out
    
    def close(self):
        self.tleap.terminate()
        del self.tleap  
    
def takeUnit(OFHandler, unitID):
    """
    This function returns the corresponding unit instance in the OFHandler giving an integer or a unit name.
    """
    try:
        unitID = int(unitID)
        return OFHandler.all.values()[unitID]
    
    except ValueError:
        return OFHandler[unitID]