import os
import sys
import time
import re
import copy
import fnmatch

import numpy as npy

try:
    import Biskit
except:
    print >>sys.stderr,"Biskit not loaded. Some functionalities may not be available."
    pass   

try:
    import Biskit.BiskitFF as bff
except:
    print >>sys.stderr,"BiskitFF not loaded. Some functionalities may not be available."
    pass



class GRID:
  """Class to manage grid files"""
  
  def __init__(self):
    pass
  
  
  def __getitem__(self, point):
    """Iterate Data giving cartesian coordinates"""
    index = self.getIndex(point)
    return self.data[index]
  
  
  def __setitem__(self, key, value):
    """When used as self[a,b,c]=x for a,b,c as cartesian coordinates
    Thus we can assign grid values giving cartesian coordinates instead of indexes."""
    index = self.getIndex(key)
    self.data[index] = value
  
  
  def __repr__(self):
    if self.source =='': return "Empty Grid instance"
    return "Grid instance for %s"%(self.source)
  
  
  def __str__(self):
    if self.source =='': return "Empty Grid instance"
    return "Grid instance for %s"%(self.source)
  
  
  def __call__(self):
    """Returns information about the grid if already defined"""
    if (self.data !='' and self.origin !='' and self.delta !=''):
      print "GRID instance source:",self.source
      print "Data shape: ",self.data.shape," number of entries: ",self.data.size
      print "Grid Origin: ",self.origin
      print "Grid Spacing: ",self.delta
    else: print "Empty grid instance or data missing!"
  
  
  def loadData(self, pick):
    """Loads grid data from existing pickle"""
    if os.path.exists(pick):
        grid = npy.load(pick)
        self.data = grid
        self.source = pick
    else: print "File not found."
  
  
  def loadVar(self, pick):
    """Loads origin and delta information from a pickled file"""
    if os.path.exists(pick):
        V = npy.load(pick)
        origin, delta= V[:]
        self.origin = origin
        self.delta = delta
    else: print "File not found."
    
  
  def writeDX(self, dxname):
    """
    Writes data into a DX Formated File
    """
    if (self.data !='' and self.origin !='' and self.delta !=''):
        grid = self.data
        origin = self.origin
        delta = self.delta
        dxf = open(dxname,'w')
        dxf.write("""#Grid generated with GRID.py\n# Time: %s\n# Source: %s\nobject 1 class gridpositions counts %i %i %i\norigin %.2f %.2f %.2f\ndelta %.5f   0   0\ndelta 0   %.5f   0\ndelta 0   0   %.5f\nobject 2 class gridconnections counts %i %i %i\nobject 3 class array type double rank 0 items %i data follows\n"""%(time.ctime(),self.source,grid.shape[0],grid.shape[1],grid.shape[2],origin[0],origin[1],origin[2],delta[0],delta[1],delta[2],grid.shape[0],grid.shape[1],grid.shape[2],grid.size))
        
        count=0
        for data in grid.flat:
            if count==3: dxf.write('\n'); count=0
            dxf.write('%6.3f\t'%(float(data)))
            count+=1
        dxf.write('\n')
        dxf.close()
        
    else: print "Data or parameters missing. Can't write DX File."
  
  
  
  def dump(self, dataname=None, varname=None):
    """This method uses by default the source name of self to generate the dumped file's names. For a specific name, it must be provided when calling the method."""
    if (self.data !='' and self.origin !='' and self.delta !=''):
      
      if not self.source:
          return "Must provide a name for the dumped files. Source attribute not defined."      
      
      if dataname:
        dname = dataname      #Grid data pickle name
      else:
        dname = str(os.path.splitext(self.source)[0])+"_grid.pick"   
      
      if varname:
        vdname = varname    #Grid variables pickle name
      else:
        vdname = str(os.path.splitext(self.source)[0])+"_grid_var.pick"
      
      npy.save(dname, self.data)
      V = [self.origin, self.delta]
      npy.save(vdname, V)
      print "Dumped pickles: ",dname," and ",vdname
    else:
      print "Data or parameters missing. Check attributes."
  
  
  
  def writeXPLOR(self, xplorname):
    """
    Write data into XPLOR formated file.
    """
    if (self.data !='' and self.origin !='' and self.delta !=''):   #If all parameters need are defined
        
        grid = self.data
        origin = self.origin
        delta = self.delta
        if isinstance(delta, npy.ndarray):
            delta = delta.tolist()
        elif type(delta) is int or type(delta) is float:
            delta = map(float, [delta, delta, delta])
        nx, ny, nz = grid.shape
        
        #minx, miny, minz, = origin[0]/delta[0], origin[1]/delta[1], origin[2]/delta[2]
        #maxx, maxy, maxz = minx+nx-1, miny+ny-1, minz+nz-1
        mind = npy.round(npy.array(origin) / delta)
        maxd = npy.round(mind + grid.shape - 1)
        space = npy.array(delta)*grid.shape   #Calculate total system expansion (distance) in x, y, z.
        deltax, deltay, deltaz = 90., 90., 90.      #Assuming rectangular shape
        
        of = open(xplorname, 'w')
        #Header and grid format
        of.write("""XPLOR Format\n       1\nFree binding energies GRID generated with count2DG.py\n%8d%8d%8d%8d%8d%8d%8d%8d%8d\n%12.3f%12.3f%12.3f%12.3f%12.3f%12.3f\nZYX\n"""%(nx, mind[0], maxd[0], ny, mind[1], maxd[1], nz, mind[2], maxd[2], space[0], space[1], space[2], deltax, deltay, deltaz))
        
        #Write Grid data
        for z in xrange(nz):
            of.write("%8d\n"%(mind[2]))
            for y in xrange(ny):
                nl = 0
                for el in grid[:,y,z]:
                    nl+=1
                    of.write("%12.5f"%(el),)
                    if nl == 6: of.write("\n"); nl=0
                if (nx%6 != 0): of.write("\n")
            mind[2]+=1
        of.close()
        
    else: print "Data or parameters missing. Can't write XPLOR File."
  
  
  def averageData(self, contract=True):
    """ 
    Will average each grid point by the surrounding ones.
    Points in the edges will be removed, thus size of the grid will change
    by one point. Also origin will be modified accordingly.
    By default, the edging points will be removed. Use contract FALSE to avoid this
    """
    
    calcgrid = npy.zeros((3,3,3))
    calcgrid += 1/27. 
    
    narr = copy.deepcopy(self.data)
    for x in xrange(1,narr.shape[0]-1):
        for y in xrange(1,narr.shape[1]-1):
            for z in xrange(1,narr.shape[2]-1):
                sl = self.data[x-1:x+2,y-1:y+2,z-1:z+2]
                narr[x,y,z] = (sl * calcgrid).sum()
    
    # Update data
    self.update(narr)
    
    if contract:
    	# Remove borders. Origin is modified in contract method already.
    	maxdelta = npy.array(self.delta).max()
    	self.contract(maxdelta)
    
  
  def count2DG(self, bvol, natoms, snaps, T, correction=1):
    
    """
    =============Convert Count Grids to Free Energy =================================
       bvol    -> Equilibrated solvent box volume
       natoms  -> Number of atoms (of specific atom type) in the initial equilibrated box
       snaps   -> Number of snapshots in the simulation
       T       -> Temperature of the simulation
       correction -> If there is a ratio solvent/water variation. This number shall correct it. 
                     By default it is 1, no correction applied.
                     
                     ---- DGBind calulation ----------
                      DGi = -kb*T*ln(Ni/No) kcal/mol
    """
    
    Cgrid = self.data   # Grid count data
    Kb = 1.987/1000.    # Boltzmann Constant    
    bvol, T = map(float, (bvol, T) )
    natoms, snaps = map(int, (natoms, snaps) )

    # print "Calculating Binding Free Energies (DGBind)..."
    
    # PROBABILITY CALCULATION
    # Given number of the atom type per molecule in the initial solvent box, 
    # and the box volume, we calculate the probability that one atom falls into a single grid point
    cell_volume = 1. / npy.prod(self.delta)
    N0 = ((natoms / bvol) / cell_volume) * snaps * correction
    
    # Convert zeros to ones to avoid log(0)
    # later will restore values after applying formula
    maskZ = Cgrid == 0   
    Cgrid[maskZ] = 1
    DGgrid = - Kb * T * ( npy.log(Cgrid) - npy.log(N0) )
    DGgrid[maskZ] = 1
    
    # print "DONE"
    return DGgrid
    
  
  def update(self, data):
    """Updates data information"""
    self.data = data
  
  
  
  def getCartesian(self, indexes):
    """Returns cartesian coordinates of the index values passed as arguments."""
    if len(indexes)==3:
      indexes = npy.array(indexes, dtype=float)
      spacing = self.delta
      origin = self.origin
      return (indexes * spacing) + origin
    else: print "Error with indexes tuple or list."
  
  
  
  def getIndex(self, cartesian):
    """Returns grid indexes for a given cartesian point (x,y,z)"""
    
    spacing = self.delta
    origin = self.origin
    if type(cartesian) is not npy.array: cartesian = npy.array(cartesian)
    ind = npy.floor((cartesian - origin) / spacing)
    # Check that index obtained is inside the grid.
    if ( npy.any( ind < 0. ) or npy.any(ind >= self.data.shape ) ):    
        return 999.
    else:
        return tuple(ind.astype(int).tolist())
 
  
  def cancelPoints(self, point, cutoff, value=999 ):
    """
    Set a grid point, and the sorrounding radio(given by cutoff in Angstroms)
    a certain value to anull it.
    Usage: self.cancelPoints([i,j,k], cutoff_in_Angstroms, value_to_assign(default=999))
    """
    data = self.data
    delta = npy.array(self.delta)
    i,j,k = point
    cutoff = npy.round(cutoff / delta).astype(int)
    
    # Avoid negative indexing when cutoff is bigger than some point index
    maskneg = point < cutoff.min()
    if npy.any(maskneg):
      cutoff[maskneg] = npy.array(point)[maskneg]
    
    #print cutoff, point, self.getCartesian(point)
    #print data[ i-cutoff[0]:i+cutoff[0]+1 , j-cutoff[1]:j+cutoff[1]+1 , k-cutoff[2]:k+cutoff[2]+1 ].min()
    data[ i-cutoff[0]:i+cutoff[0]+1 , j-cutoff[1]:j+cutoff[1]+1 , k-cutoff[2]:k+cutoff[2]+1 ] = value
    #print data[ i-cutoff[0]:i+cutoff[0]+1 , j-cutoff[1]:j+cutoff[1]+1 , k-cutoff[2]:k+cutoff[2]+1 ].min()
    
  
  
  
  def mergeDeleteProt(self, PDB, cutoff, value=999):
    
    """
    Method to "Delete" grid information on the points overlaping atoms of a PDB.
    
    @param1 PDB:     protein
    @param1 type:    Biskit PDBModel instance or PDB File name
    
    @param2 cutoff: Cutoff allows to expand the assignement to N angstroms around that point.
    @param2 type;   int or float
    
    @param3 value:  Optionally one can provide the value to assign to those points (default = 999)
    @param3 type:   int
    """
    
    if isinstance(PDB, Biskit.PDBModel):
        protein = PDB
        
    else:
        protein = Biskit.PDBModel(PDB)
        
    origin = self.origin
    delta = self.delta
    
    for atom in protein.xyz:
      idx = self.getIndex(atom)
      #a,b,c = atom
      #i, j, k = int(round((a-origin[0])/delta[0])), int(round((b-origin[1])/delta[1])), int(round((c-origin[2])/delta[2]))
      if idx != 999.: self.cancelPoints(idx, cutoff, value)
  
  
  
  def mergeConserveProt(self, PDB, cutoff, value=999):
    
    """
    Method to Conserve only grid information near the points overlaping atom's coordinates of a PDB.
    
    @param1 PDB:     protein
    @param1 type:    Biskit PDBModel instance or PDB File name
    
    @param2 cutoff: How far from the atom coordinate should keep grid information. In Angstroms.
    @param2 type;   int or float
    
    @param3 value:  Value to assign to the rejected points (default = 999)
    @param3 type:   int
    
    return data
    """
        
    if isinstance(PDB, Biskit.PDBModel):
        protein = PDB
        
    else:
        protein = Biskit.PDBModel(PDB)
        
    origin = self.origin
    delta = self.delta
    _dataBackup = copy.deepcopy(self.data)
    self.data = npy.zeros(self.data.shape) == 0
    
    for atom in protein.xyz:
      a,b,c = atom
      i, j, k = int(round((a-origin[0])/delta[0])), int(round((b-origin[1])/delta[1])), int(round((c-origin[2])/delta[2]))
      self.cancelPoints([i,j,k], cutoff, False)
    
    _dataBackup[self.data] = value
    
    return _dataBackup
  
  
  def copy(self):
    """Returns a new instance as a deepcopy of self"""
    new_inst = copy.deepcopy(self)
    return new_inst
    
  
  def min_indx(self):
    """Returns index for the minimum value in data"""
    return npy.array( zip( *npy.where( self.data == self.data.min() ) ) )[0]
  
  
  def expand(self, buff,fill=0):
    """
    Expand the array size building a Zero filled Array adding buff Angstrom each side.
    Place data inside and modify origin information.
    """
    delta = self.delta
    size = self.data.shape
    buffering = [ int( buff / delta[0] ), int( buff / delta[1] ), int( buff / delta[2] ) ]
    
    G = npy.zeros( [ size[0] + (buffering[0]*2), size[1] + (buffering[1]*2), size[2] + (buffering[2]*2) ] )+fill
    G[ buffering[0]:size[0] + buffering[0], buffering[1]:size[1] + buffering[1], buffering[2]:size[2] + buffering[2] ] = self.data.copy()
    
    self.data = G.copy()
    self.origin = npy.array(self.origin) - buff
  
    
  def contract(self, buff):
    """ 
    Inverser operation to expand. 
    buff in angstroms.
    """
    delta = self.delta
    size = self.data.shape
    buffering = [ int( buff / delta[0] ), int( buff / delta[1] ), int( buff / delta[2] ) ]
    self.data = self.data[ buffering[0]:size[0]-buffering[0], buffering[1]:size[1]-buffering[1], buffering[2]:size[2]-buffering[2] ].copy()
    self.origin = npy.array(self.origin) + buff
  
  
  def trim(self, X):
    """
    Trim self array data to fit common data with grid X. 
    New shape and coordinates is common in the two arrays but data is copied from self
    """
    if self.delta == X.delta:
        
        A = self
        
        Amaxind = npy.array( A.data.shape )     
        Xmaxind = npy.array( X.data.shape )
        
        #Get new Origin and new Shape
        newOrigin = npy.vstack((A.origin, X.origin)).max(axis=0)
        newMaximum = npy.vstack((A.getCartesian(Amaxind), X.getCartesian(Xmaxind))).min(axis=0)
        newShape = ( newMaximum - newOrigin ) / A.delta
        newShape.astype('int')
        
        newarray = create(newShape, newOrigin, A.delta)
        newarray.source = 'Trimmed array'
        
        #Copy self data to the new array
        init = A.getIndex(newarray.origin)
        newarray.update( A.data[ init[0]:init[0]+newarray.data.shape[0], init[1]:init[1]+newarray.data.shape[1], init[2]:init[2]+newarray.data.shape[2] ].copy() )
        
        return newarray
            
    else:
        print "Spacing is not equal. Can not trim."

  
  def takeSubGridBox(self, bot, top,forceCubic=0):
    
    """
    Return a grid instance containing a part of the main grid.
    
    arg1:   bot :   bottom Cartesian coordinates of the subgrid
    arg1    type:   tuple or list [x,y,z]
    
    arg2:   top:    top Cartesian coordinates of the subgrid
    arg2    type:   tuple or list [x,y,z]

    arg3:   forceCubic :    a boolean to know if we force  cubic grid to be returned
    arg3:   type:   int (0 or 1)
    """

    subOrigin = bot
    subOriginIndx = npy.array(self.getIndex(bot))
    subTopIndx = npy.array(self.getIndex(top))
    
    if forceCubic:
        diffIndx=subTopIndx-subOriginIndx
        maxDiff=npy.max(diffIndx)
        subTopIndx+=maxDiff-diffIndx

    subGrid = GRID()
    subGrid.delta = self.delta
    subGrid.origin = subOrigin
    tmp=self.data[ subOriginIndx[0]:subTopIndx[0]+1, subOriginIndx[1]:subTopIndx[1]+1, subOriginIndx[2]:subTopIndx[2]+1 ]
    subGrid.data = npy.ones_like(tmp)
    subGrid.data[:] = tmp[:]
    
    subGrid.source = self.source + 'SUBGRID'
    
    return subGrid
    
  
  
  def takeSubGridPoint(self, point, distance):
    
    """
    Return a grid instance containing a part of the main grid.
    arg1:   point:    coordinates for the center of the box
    arg1    type:      tuple or list [x,y,z]
    
    arg2:   distance:   Distance of the box side size in Angstroms
    arg2    type:       tuple or list [dx, dy, dz] for rectangular box
    arg2    type:       int or float for cubic box
    """
    
    subOrigin = npy.array(point) - (distance/2.)
    subOriginIndx = self.getIndex(tuple(subOrigin))
    subTop = npy.array(point) + (distance/2.)
    subTopIndx = self.getIndex(tuple(subTop))
    
    subGrid = GRID()
    subGrid.delta = self.delta
    subGrid.origin = subOrigin
    subGrid.data = self.data[ subOriginIndx[0]:subTopIndx[0]+1, subOriginIndx[1]:subTopIndx[1]+1, subOriginIndx[2]:subTopIndx[2]+1 ]
    subGrid.source = self.source + 'SUBGRID'
    
    return subGrid
  
  def getSphereValues(self, center, r):
        """ This method returns an array containing the values around a certain
        cartesian point given a radius """
        
        idx = self.getIndex(center)
        spacing = npy.array(self.delta).max()
        span = npy.round(float(r) / spacing)
        
        return self.data[idx[0]-span:idx[0]+span+1, idx[1]-span:idx[1]+span+1, idx[2]-span:idx[2]+span+1]


class create_DEPRECATED(GRID):
    def __init__(self, gridShape, origin, spacing, value_filling=0):
        """
        Create a empty or value_filled grid with all the given parameters. Dimensions in Angstroms in x,y,z.
        """
        print gridShape
        print spacing
        print origin
        if (len(gridShape)==3 and len(spacing)==3 and len(origin)==3):
            self.data = npy.zeros(gridShape) + value_filling
            self.origin = origin
            self.delta = spacing
            self.source = "User created GRID"
            self.shape = [int(el) for el in gridShape]
        else: print "Bad arguments. Cannot create grid."
    

class create(GRID):
    def __init__(self, gridShape, origin, spacing=None, geomExtent=None, value_filling=0, margin=0.0, **kwargs):
        """
        Create a empty or value_filled grid with all the given parameters. Dimensions in Angstroms in x,y,z.
        """
        #print gridShape
        #print spacing
        #print origin
        
        if spacing is not None:
            self.data = npy.zeros(gridShape, **kwargs) + value_filling
            self.origin = origin
            self.source = "User created GRID"
            self.shape = map(int, gridShape)
            if type(spacing) is int or type(spacing) is float:
                self.delta = [spacing, spacing, spacing]
            else:
                self.delta = spacing
        else :
            print margin
            self.origin=origin-margin
            self.delta=(geomExtent+margin*2)/gridShape
            self.shape= map(int, gridShape)
            self.data = npy.zeros(gridShape, **kwargs) + value_filling
            self.source = "User created GRID"
            

class createFromPDB(GRID, create): 
    def __init__(self, PDB, spacing, buff=0, value_filling=0, takeProtein=True):
	"""
	PDB	PDBModel or PDB file
	spacing grid spacing
	buff	buffer around atom coordinates in angstroms
	value_filling	initial value of the grid
	takeProtein	Boolean. If true, will try to compress only protein from PDB.
			If False, will take the whole PDB itself
	"""
        if isinstance(PDB, Biskit.PDBModel):
            prot = PDB
        elif os.path.exists(PDB):
            pdb_model = Biskit.PDBModel(PDB)                    #Open pdb as Biskit PDBModel
            if takeProtein: prot = pdb_model.compress(pdb_model.maskProtein())  #Read only protein from PDBModel
	    else: prot = pdb_model
        else: 
            print "Couldn't locate PDB file: ",PDB
            return False        
        
        xyz = prot.xyz                                                  # get protein coordinates
        extension = ( xyz.max(axis=0) - xyz.min(axis=0) ) + buff * 2    # protein extension in x,y,z axis
        shape = (extension / spacing).astype(int)                       # New grid shape
        origin = xyz.min(axis=0) - buff               # grid min cartesian coordinates as origin 
        
        create.__init__(self, shape, origin.tolist(), spacing=spacing, value_filling=value_filling)
        self.source = 'Generated from PDB: '+str(PDB)


class createFromMol2(create): 
  def __init__(self, Mol2Model, shape, margin=0, value_filling=0):
        
        xyz = Mol2Model.xyz                                             #get molecule coordinates
        extension = ( xyz.max(axis=0)-xyz.min(axis=0) ) + margin*2      #protein extension in x,y,z axis
        spacing = extension/npy.array(shape)                             #New grid size with buff*2 angstroms added in each axis
        origin = xyz.min(axis=0) - margin            #grid min cartesian coordinates as origin build a zeros grid with 0.5 spacing
        #print "Applying ",buff," Angstroms of buffer each side."
        #print "Origin: ",origin
        #print "Spacing: ",spacing
        create.__init__(self,shape,origin,spacing,extension, value_filling)
        self.source = 'Generated from Mol2 file: '+str(Mol2Model)

class get_grid(GRID):
  
  def __init__(self, FILE = None, buff = None):
    self.data = ''
    self.origin = ''
    self.delta = ''
    self.atype = ''
    self.source = ''
    
    if FILE:
      if fnmatch.fnmatch(FILE, '*.xplor'):
        self.parseXPLOR(FILE, buff)
        
      elif fnmatch.fnmatch(FILE, '*.dx'):
        self.parseDX(FILE, buff)
        
      else: print "Can not identify ",FILE,"file format.It should contain 'xplor' or 'dx' extension"
        
  def parseXPLOR(self, XPLOR, buff = None):
    if os.path.exists(XPLOR):
        def readBlock(fh, size, dType):
            res = npy.fromstring(fh.read(size),sep=" ",dtype=dType)
            return res.copy()
        #Open XPLOR Grid File
        f = open(XPLOR, 'r')
        
        #Skip 1st line
        f.readline()
        #Header lines:
        head = f.readline().split()
        head = int(head[0])
        while head:
            f.readline()
            head-=1
        
        #Read grid format
        format1 = f.readline().split()
        for i in xrange(len(format1)): format1[i] = int(format1[i])
        nx, minx, maxx, ny, miny, maxy, nz, minz, maxz = format1
        
        format2 = f.readline().split()
        for i in xrange(len(format2)): format2[i] = float(format2[i])
        dx, dy, dz, deltax, deltay, deltaz = format2
        
        #Skip ZYX line
        f.readline()
        
        #Calculate size to read
        #XPLOR Format has ZYX format, X fast, Y medium, Z slow.
        #X data is sorted in 6 columns
        #Each data is 12 bits and each line is 1 bit more (newline charater)
        #This format is repeated ny times. Then a new Z block starts (one line to be skipped then).
        size = int((nx*12)+npy.ceil(nx/6.))
        
        grid = npy.zeros([nx,ny,nz])
        
        for z in xrange(nz):
            f.readline() #skip block identifier line
            for y in xrange(ny):
                grid[:,y,z]=readBlock(f, size, 'float')
        
        #Prepare standard descriptors:
        delta = [(dx/nx), (dy/ny), (dz/nz)]
        origin = [minx*delta[0],miny*delta[1],minz*delta[2]]
        
        #Assign attributes
        self.data = grid
        self.origin = origin
        self.delta = delta
        self.source = XPLOR
        
        #Expand if buffer is assigned
        if buff is True:
          self.expand(buff)

    else: print "Could'nt find ",XPLOR,"file in current folder. Check path and name."

  def parseDX(self, DX, buff = None):
    """ 
    Read data from a DX file
    parmam  buff:   buffer to expand the initial grid with zeros to a given number of ansgtroms each side
    type            float
    """
    
    if os.path.exists(DX):
        f=open(DX,"r")
        #read the header - here is an example
        header=""
        for i in range(3) : header= header + f.readline()
        
        #read the grid size
        r=re.compile('\w+')
        gsize=r.findall(f.readline())
#        print gsize
        gsize=[int(gsize[-3]),int(gsize[-2]),int(gsize[-1])]
        
        #read the origin of the system
        line=f.readline().split()
        origin=[float(line[-3]),float(line[-2]),float(line[-1])]
        
        #read grid space
        line=f.readline().split()
        deltax=[float(line[-3]),float(line[-2]),float(line[-1])]
        line=f.readline().split()
        deltay=[float(line[-3]),float(line[-2]),float(line[-1])]
        line=f.readline().split()
        deltaz=[float(line[-3]),float(line[-2]),float(line[-1])]
        
        #pay attention here, this assumes always orthogonal normalized space, but normally it should be ok
        delta=[deltax[0],deltay[1],deltaz[2]]
        
        #read the number of data
        f.readline()
        r=re.compile('\d+')
        n_entries=int(r.findall(f.readline())[2])
                
        #check correpondence with expected data points
        if(n_entries!=gsize[0]*gsize[1]*gsize[2]) : sys.exit("Error reading the file. The number of expected data points does not correspond to the number of labeled data points in the header.")
        
        #load data into numpy array
        Grid = npy.fromstring(f.read(), sep=' ', dtype=float).reshape(gsize) #reshaping to fit grid format (it keeps Z fast, Y medium, X slow data organization)
        
        if Grid.size != n_entries: sys.exit("Error reading the file. The number of expected data points does not correspond to the number of labeled data points in the header.")
        f.close()
        
        self.data= Grid
        self.origin = origin
        self.delta = delta
        self.source = DX
        
        if buff != None:
          self.expand(buff)
        
    else: print "Could'nt find ",DX,"file in current folder. Check path and name."


def calculateEnergyGrid(proteinPdbFile, atomtypePdbFile, proteinTopologyFile, atomtypeTopologyFile, clearance = 2., spacing= [0.5, 0.5, 0.5]):
    prot = Biskit.PDBModel(proteinPdbFile)
    prot = prot.compress( prot.maskProtein() )
    prot.loadAmberTopology(proteinTopologyFile)
    probe = Bsikit.PDBModel(atomtypePdbFile)
    probe.loadAmberTopology(atomtypeTopologyFile)
    
    t0 = time.time()
    print "All loaded. Calculating grid..."
    grid = createFromPDB(proteinPdbFile, spacing, clearance)
    print grid()
    
    pcenter = probe.xyz.mean(axis=0)
    print pcenter
    
    for i in xrange(grid.data.shape[0]):
      for j in xrange(grid.data.shape[1]):
        for k in xrange(grid.data.shape[2]):
            print i,j,k
            gridCartesian = grid.getCartesian([i,j,k])
            translationMod = npy.array(gridCartesian) - pcenter
            probe.xyz += translationMod
            grid.data[i,j,k] = bff.getInteractionEnergyCpu(prot, probe, dict=0)[-1]
    
    t1=time.time()
    print "DONE in %.2f s"%((t1-t0)/60.)
    return grid


def tanimotoIdentityIsoval(Grid1, Grid2, isovalue, up=True, expand = False, verbose = False):
    """
    This functions returns a Tanimoto coefficient to compare two grids ressemblance.
        Tanimoto = Nab / (Na + Nb + Nab)
        Nab = Ones in A and B
        Na = Ones in A
        Nb = Ones in B
    
    @arg1   Grid1:      Grid instance
    @arg2   Grid2:      Grid instance
    @arg3   isovalue:   cut value to transform to binary data
            type:       int or float
    @arg4   up:         if up is True, take grid values over isovalue, else take values under isovalue
            type:       bool
    @arg5   expand:     Expand the 1s this number of angstroms around. This way we avoid grid spacing problem.
            type:       float or int
    @arg6   verbose:    If True, it returns also Na, Nb, Nab, minimum value and cutoff value for each grid as dictionary
            type:       bool
    """
    
    # Grids are first trimmed to obtain the common part in space.    
    A, B = trim(Grid1, Grid2)
    A.source = Grid1.source
    B.source = Grid2.source
    
    # Mask TRUE (1) values under/over isovalue
    
    if up is True:
        AmaskTrue = A.data >= isovalue
        BmaskTrue = B.data >= isovalue
    else:
        AmaskTrue = A.data <= isovalue
        BmaskTrue = B.data <= isovalue   
    
    # EXPAND POSITIVE POINTS 0.5 ANGSTROMS AROUND
    if expand is True:
        ExpandA = A.copy()
        ExpandB = B.copy()
        ExpandA.data = AmaskTrue*1
        ExpandA.expand(2)   #Expand to avoid bad indexing
        ExpandB.data = BmaskTrue*1
        ExpandB.expand(2)
        for grid in (ExpandA, ExpandB):
            point_lst = npy.vstack(npy.nonzero(grid.data==1)).T
            for point in point_lst:
                grid.cancelPoints(point = point, cutoff = expand, value = 1)
        ExpandA.contract(2)
        ExpandB.contract(2)
        AmaskTrue = ExpandA.data == 1
        BmaskTrue = ExpandB.data == 1
    
    # Mix bitwise AND
    ABTrue = AmaskTrue * BmaskTrue
    
    # Tanimoto calculation
    Na, Nb, Nab = len(npy.nonzero(AmaskTrue)[0]), len(npy.nonzero(BmaskTrue)[0]), len(npy.nonzero(ABTrue)[0])
    Tc = float(Nab) / (Na + Nb - Nab)
    
    if verbose is False:
        return Tc
    else:
        return {'TanimotoCoeff':Tc, 'Na':Na, 'Nb':Nb, 'Nab':Nab, 'Amin': A.data.min(),'Bmin':B.data.min()}


def tanimotoIdentityCutoff(Grid1, Grid2, cutoff, up=False, expand = False, verbose = False):
    """
    This functions returns a Tanimoto coefficient to compare two grids ressemblance.
        Tanimoto = Nab / (Na + Nb + Nab)
        Nab = Ones in A and B
        Na = Ones in A
        Nb = Ones in B
    
    @arg1   Grid1:      Grid instance
    @arg2   Grid2:      Grid instance
    @arg3   isovalue:   cut value to transform to binary data
            type:       int or float
    @arg4   up:         if up is True, take grid values over isovalue, else take values under isovalue
            type:       bool
    @arg5   expand:     Expand the 1s 0.5 angstroms around. This way we avoid grid spacing problem.
            type:       bool
    @arg6   verbose:    If True, it returns also Na, Nb, Nab, minimum value and cutoff value for each grid as dictionary
            type:       bool
    """
    
    #Grids are first trimmed to obtain the common part in space.    
    A, B = trim(Grid1, Grid2)
    A.source = Grid1.source
    B.source = Grid2.source
    
    Amin = A.data.min()        
    Bmin = B.data.min()   
    
    #Flatten negative data, sort it and get the value of the CUT position of negative point as cutoff
    #Sort number of points for each grid and choose the lower one
    
    fsortA = npy.array( A.data[A.data<0].flat )
    fsortA.sort()
    fsortB = npy.array( B.data[B.data<0].flat )
    fsortB.sort()    
    
    cutIndex = []
    [ cutIndex.append( int( (sortedlist.shape[0]-1) * (cutoff/100.) )  ) for sortedlist in (fsortA, fsortB) ]
    
    cutIndex.sort()
    finalIndex = cutIndex[0]
    
    AcutValue = fsortA[finalIndex]
    BcutValue = fsortB[finalIndex]    
    
    #Mask TRUE (1) values under/over isovalue
    
    if up is True:
        AmaskTrue = A.data > AcutValue
        BmaskTrue = B.data > BcutValue
    
    else:
        AmaskTrue = A.data < AcutValue
        BmaskTrue = B.data < BcutValue  
    
    #EXPAND POSITIVE POINTS 0.5 ANGSTROMS AROUND
    if expand is True:
        ExpandA = A.copy()
        ExpandB = B.copy()
        ExpandA.data = AmaskTrue*1
        ExpandA.expand(2)   #Expand to avoid bad indexing
        ExpandB.data = BmaskTrue*1
        ExpandB.expand(2)
        for grid in (ExpandA, ExpandB):
            point_lst = npy.vstack(npy.nonzero(grid.data==1)).T
            for point in point_lst:
                grid.cancelPoints(point = point, cutoff = 0.5, value = 1)
        ExpandA.contract(2)
        ExpandB.contract(2)
        AmaskTrue = ExpandA.data == 1
        BmaskTrue = ExpandB.data == 1
    
    Na, Nb = AmaskTrue.sum(), BmaskTrue.sum()
    
    #Mix bitwise AND
    Nab = (AmaskTrue * BmaskTrue) *1
    Nab = Nab.sum()

    Tc = float(Nab) / (Na + Nb - Nab)
    
    if verbose is False:
        return Tc
    else:
        return {'TanimotoCoeff':Tc, 'Na':Na, 'Nb':Nb, 'Nab':Nab, 'Amin': Amin,'Bmin':Bmin}



def similarityIndexIsoval(Grid1, Grid2, isovalue, up = True, expansion = 1., verbose = False):
    """
    This functions returns a modified Tanimoto coefficient to compare two grids ressemblance.
    It takes in account the expansion in the calculation.
        a' = expanded a without original points. Thus only expansion.
        b' = expanded b without 1s in b
        SimIndex = (Nab + Na'b + Nab' + Na'b') / (Na + Nb + Na'b' + (Na'b/2) + (Nab'/2) - Nab)
        Na = Ones in original A
        Na'= Ones in expansion of A
        Nab = Ones common in A and B        
        Na'b = Ones in common in expansion of A and original B
        etc..
    
    @arg1   Grid1:      Grid instance
    @arg2   Grid2:      Grid instance
    @arg3   isovalue:   cut value to transform to binary data
            type:       int or float
    @arg4   up:         if up is True, take grid values over isovalue, else take values under isovalue
            type:       bool
    @arg5   expansion:  Amount of expansion of true points in Angstroms.This way we avoid grid spacing problem.
            type:       int or float (default = 1)
    @arg6   verbose:    If True, it returns a dictionary with index details.
            type:       bool
            
    """
    
    #Grids are first trimmed to obtain the common part in space.    
    A, B = trim(Grid1, Grid2)
    A.source = Grid1.source
    B.source = Grid2.source
    
    Amin = A.data.min()
    Bmin = B.data.min()
    
    #Mask TRUE (1) values under/over isovalue
    
    if up is True:
        Amask = A.data >= isovalue
        Bmask = B.data >= isovalue
    
    else:
        Amask = A.data <= isovalue
        Bmask = B.data <= isovalue
    
    #EXPAND POSITIVE POINTS 0.5 ANGSTROMS AROUND
    A_ = A.copy()
    B_ = B.copy()
    A_.data = Amask*1
    A_.expand(2)   #Expand to avoid bad indexing
    B_.data = Bmask*1
    B_.expand(2)
    for grid in (A_, B_):
        point_lst = npy.vstack(npy.nonzero(grid.data==1)).T
        for point in point_lst:
            grid.cancelPoints(point = point, cutoff = expansion, value = 1)
    A_.contract(2)
    B_.contract(2)
    
    A = Amask*1
    B = Bmask*1
    
    A_Mask = (A_.data - A)==1
    B_Mask = (B_.data - B)==1
    
    #Mix bitwise AND and sum to obtain number of 1s
    A  = float( A.sum() )
    B  = float( B.sum() )
    A_ = float( A_.data.sum() )
    B_ = float( B_.data.sum() )
    AB  = float( ( ( Amask*Bmask)*1  ).sum() )
    A_B = float( ( (A_Mask*Bmask)*1  ).sum() )
    AB_ = float( ( (Amask*B_Mask)*1  ).sum() )
    A_B_= float( ( (A_Mask*B_Mask)*1 ).sum() )
    
    #SimilarityIndex calculation
    SimIdx = ( AB + A_B + AB_ + A_B_) / (A + B + A_B_ + (A_B/2.) + (AB_/2.) - AB)
    
    if verbose is False:
        return SimIdx
   
    else:
        return {'SimilarityIndex':SimIdx, 'A':A, 'B':B, 'AB':AB, 'A_':A_, 'B_':B_, 'A_B':A_B, 'AB_':AB_, 'A_B_':A_B_, 'Amin': Amin, 'Bmin':Bmin}


def _cutValue(A, B, cutoff):
    """
    Complement function to trim grids.
    It returns the value for A and B grids that will serve as 
    cutoff % threshold to return approximately that amount of points
    """

    #Flatten negative data, sort it and get the value of the CUT position of negative point as cutoff
    #Sort number of points for each grid and choose the lower one
    
    fsortA = npy.array( A.data[A.data<0].flat )
    fsortA.sort()
    fsortB = npy.array( B.data[B.data<0].flat )
    fsortB.sort()    
    
    cutIndex = []
    [ cutIndex.append( int( (sortedlist.shape[0]-1) * (cutoff/100.) )  ) for sortedlist in (fsortA, fsortB) ]
    
    cutIndex.sort()
    finalIndex = cutIndex[0]
    
    AcutValue = fsortA[finalIndex]
    BcutValue = fsortB[finalIndex]    

    return AcutValue, BcutValue


def similarityIndex2(Grid1, Grid2, cutoff= 5., verbose = False, trimfirst=True):
    """
    Obtain Similarity Comparison using SimIndex.
    Both grids will be trimmed to contain same number of points.
    
    cutoff (float)  -   Threshold of most negative points to be compares.
                        E.g. 5% of most negative points.
                        It will take into account the number of points independently
                        of the threshold value.
        
    trimfirst (bool)    If True, cutting value will be obtained before trimming.
                        If False, cutting value will be obtained after trimming (thus
                        we only consider values in the trimmed zone).
    
    """
    if not trimfirst: AcutValue, BcutValue = _cutValue(Grid1, Grid2, cutoff)
        
    A, B = trim(Grid1, Grid2)
    A.source = Grid1.source
    B.source = Grid2.source
    
    Amin = A.data.min()        
    Bmin = B.data.min()   
    
    if trimfirst: AcutValue, BcutValue = _cutValue(A, B, cutoff)    
    
    #Mask TRUE (1) values under cutoff
    Amask = A.data < AcutValue     
    Bmask = B.data < BcutValue
    
    #EXPAND POSITIVE POINTS 0.5 ANGSTROMS AROUND
    A_ = A.copy()
    B_ = B.copy()
    A_.data = Amask*1
    A_.expand(2)   #Expand to avoid bad indexing
    B_.data = Bmask*1
    B_.expand(2)
    for grid in (A_, B_):
        point_lst = npy.vstack(npy.nonzero(grid.data==1)).T
        for point in point_lst:
            grid.cancelPoints(point = point, cutoff = 0.5, value = 1)
    A_.contract(2)
    B_.contract(2)
    
    A = Amask*1
    B = Bmask*1
    
    A_Mask = (A_.data - A)==1
    B_Mask = (B_.data - B)==1
    
    
    #Mix bitwise AND and sum to obtain number of 1s
    A  = float( A.sum() )
    B  = float( B.sum() )
    A_ = float( A_.data.sum() )
    B_ = float( B_.data.sum() )
    
    AB  = float( ( ( Amask*Bmask)*1  ).sum() )
    A_B = float( ( (A_Mask*Bmask)*1  ).sum() )
    AB_ = float( ( (Amask*B_Mask)*1  ).sum() )
    A_B_= float( ( (A_Mask*B_Mask)*1 ).sum() )
    
    #SimilarityIndex calculation
    SimIdx = ( AB + A_B + AB_ + A_B_) / (A + B + A_B_ + (A_B/2.) + (AB_/2.) - AB)
    
    if verbose is False:
        return SimIdx
    
    else:
        return {'SimilarityIndex':SimIdx, 'A':A, 'B':B, 'AB':AB, 'A_':A_, 'B_':B_, 'A_B':A_B, 'AB_':AB_, 'A_B_':A_B_, 'Amin': Amin, 'Acut':AcutValue, 'Bmin':Bmin, 'Bcut':BcutValue}


def trim(*Glist):
    """
    Obtain common part of all grids given as argument. 
    arg1:   Glist   List containing more than 1 GRID instances to crop.
    arg1:   type    List or several GRID instances
    
    returns a list of trimmed grid instances
    """     
    if len(Glist)<2:
        if type(Glist[0]) is list:
          Glist = Glist[0]
        else:
          return "ERROR. Must provide at least 2 grids."
        
    deltaList = npy.array( [grid.delta for grid in Glist] )
    
    if npy.any(deltaList != deltaList[0]): 
      return "ERROR. All grids should have same spacing."
    
    delta = deltaList[0]
    originList = npy.array( [grid.origin for grid in Glist] )
    maxCoordList = npy.array( [grid.getCartesian(grid.data.shape) for grid in Glist] )
    
    #Get maximum Origin and Cartesian Coordinates for minimum shape
    newOrigin = originList.max(axis=0)
    newMaximum = maxCoordList.min(axis=0)
    newShape = ( newMaximum - newOrigin ) / delta
    newShape.astype('int')
    
    trimmedList = []
    
    for grid in Glist:
        newarray = create(newShape, newOrigin, delta)
        newarray.source = 'Trimmed array of %s'%(grid.source)
        #Copy data to the new array
        init = grid.getIndex(newarray.origin)
        newarray.update( grid.data[ init[0]:init[0]+newarray.data.shape[0], init[1]:init[1]+newarray.data.shape[1], init[2]:init[2]+newarray.data.shape[2] ].copy() )
        trimmedList.append(newarray)
        
    return trimmedList


import Biskit.tools as tools

if __name__ == '__main__':
  proteinPdb = tools.testRoot()+os.sep+"Grid"+os.sep+"protein2.pdb"
  proteinTop = tools.testRoot()+os.sep+"Grid"+os.sep+"protein2.top"
  probePdb = tools.testRoot()+os.sep+"Grid"+os.sep+"methane.pdb"
  probeTop = tools.testRoot()+os.sep+"Grid"+os.sep+"methane.top"
  grid = calculateEnergyGrid( proteinPdb, probePdb, proteinTop,  probeTop )
  grid.writeDX(tools.testRoot()+os.sep+"Grid"+os.sep+"Methane_test_grid.dx")
  print grid()
  print "min",grid.data.min()
  print "max",grid.data.max()
  grid.dump('Methane_test_grid')
  
