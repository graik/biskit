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
general purpose utility methods.
"""

import sys, string
import os.path as osp
import shutil
import operator
import os, pickle  ## for Load and Dump
import tempfile
import traceback
from inspect import getframeinfo
from time import localtime
import glob
import subprocess
import gzip
import collections
import types
from biskit.core import oldnumeric as N0

class ToolsError( Exception ):
    pass

class PickleError( ToolsError ):
    pass

def errWriteln(s):
    """
    print s to standard error with line feed.
    
    :param s: string
    :type  s: str
    """
    sys.stderr.write(s+'\n')
    sys.stderr.flush()


def errWrite(s):
    """
    print s to standard error.

    :param s: string
    :type  s: str    
    """
    sys.stderr.write(s)
    sys.stderr.flush()


def flushPrint(s):
    """
    print s without line break and flush standard out.

    :param s: string
    :type  s: str    
    """
    sys.stdout.write(s)
    sys.stdout.flush()


def lastError():
    """
    Collect type and line of last exception.
    
    :return: '<ExceptionType> in line <lineNumber>:<Exception arguments>'
    :rtype: String
    """
    try:
        trace = sys.exc_info()[2]
        why = sys.exc_info()[1]
        try:
            why = sys.exc_info()[1].args
        except:
            pass
        file = getframeinfo( trace.tb_frame )[0]

        result = "%s in %s line %i:\n\t%s." % ( str(sys.exc_info()[0]),
                  file, trace.tb_lineno, str(why) )

    finally:
        trace = None

    return result


def lastErrorTrace( limit=None ):
    tb = sys.exc_info()[2]

    lines = traceback.extract_tb( tb, None )

    result = ''
    for l in lines:
        pyFile = stripFilename( l[0] )
        result += '%s: %i (%s) %s\n' % (pyFile, l[1],l[2],l[3])

    return result


def dictAdd( dic, key, value, forceList=False ):
    """
    Add value to dic, create list, if dic has already value in key.

    :param key: dictionary key
    :type  key: str
    :param value: value
    :type  value: any
    """
    if key in dic:
        old = dic[key]

        if type( old ) != list and value != old:
            dic[ key ] = [ old ] + [ value ]
        else:
            if type( old ) == list and value not in old:
                dic[ key ] = old + [ value ]

    else:
        if forceList and (type( value ) != list):
            dic[key] = [ value ]
        else:
            dic[key] = value


def absfile( filename, resolveLinks=1 ):
    """
    Get absolute file path::
      - expand ~ to user home, change
      - expand ../../ to absolute path
      - resolve links
      - add working directory to unbound files ('ab.txt'->'/home/raik/ab.txt')

    :param filename: name of file
    :type  filename: str
    :param resolveLinks: eliminate any symbolic links (default: 1)
    :type  resolveLinks: 1|0
    
    :return: absolute path or filename
    :rtype: string

    :raise ToolsError: if a ~user part does not translate to an existing path
    """
    if not filename:
        return filename
    r = osp.abspath( osp.expanduser( filename ) )

    if '~' in r:
        raise ToolsError('Could not expand user home in %s' % filename)

    if resolveLinks:
        r = osp.realpath( r )
    r = osp.normpath(r)
    return r


def homefile( filename, otherUser=1, ownCopy=1 ):
    """
    Relativize a file name to ~ or, if it is in another user's home,
    to ~otheruser or, if it is in nobody's home, to / .
    
    :class:`splithome()` is used to also guess home directories of other users.
    
    :param filename: name of file
    :type  filename: str
    :param otherUser: look also in other user's home directories (default 1)
    :type  otherUser: 1|0
    :param ownCopy: replace alien path by path into own home directory if
                    possible, e.g. ~other/data/x is replaced
                    by ~/data/x if there is such a file. (default 1) Careful!
    :type  ownCopy: 1|0
                
    :return: path or filename
    :rtype: str
    """
    f = absfile( filename )
    my_home = osp.expanduser('~')
    user_home, rest = splithome( f )

    if user_home == my_home:
        return f.replace( user_home, '~', 1 )

    if otherUser and user_home != '':
        ## first try to find same path in own home directory
        if ownCopy:
            my_path = os.path.join( my_home, rest )
            if osp.exists( my_path ):
                return my_path

        user = osp.split( user_home )[-1]
        return f.replace( user_home+'/', '~' + user + '/', 1 )

    return f


def splithome( filename ):
    """
    Split path into home directory and remaining path. Valid home directories
    are folders belonging to the same folder as the current user's home. I.e.
    the method tries also to guess home directories of other users.

    :param filename: name of file
    :type  filename: str
    
    :return: home folder of some user, remaining path relative to home
    :rtype: (str, str)
    """
    home = osp.expanduser( '~' )
    home_base = osp.split( home )[0]

    if filename.find( home_base ) != 0:
        return '', filename

    f = filename.replace( home_base + '/', '', 1 )
    user = f.split( '/' )[0]

    user_home = os.path.join( home_base, user )
    rest = f.replace( user + '/', '', 1 )

    return user_home, rest


def __pathsplit(p, rest=[]):
    """
    from ASPN Python Cookbook
    """
    (h,t) = os.path.split(p)
    if len(h) < 1: return [t]+rest
    if len(t) < 1: return [h]+rest
    return __pathsplit(h,[t]+rest)

def __commonpath(l1, l2, common=[]):
    """
    from ASPN Python Cookbook
    """
    if len(l1) < 1: return (common, l1, l2)
    if len(l2) < 1: return (common, l1, l2)
    if l1[0] != l2[0]: return (common, l1, l2)
    return __commonpath(l1[1:], l2[1:], common+[l1[0]])

def relpath(p1, p2):
    """
    Translate p2 into a path relative to p1.
    
    :param p1: base path
    :type p1: str
    :param p2: target path
    :type p2: str

    :return: relative path p1 -> p2
    :rtype: str
    """
    if not p1 or not p2:
        return p2
    (__common,l1,l2) = __commonpath(__pathsplit(p1), __pathsplit(p2))
    p = []
    if len(l1) > 0:
        p = [ '../' * len(l1) ]
    p = p + l2
    return os.path.join( *p )


def stripSuffix( filename ):
    """
    Return file name without ending.

    :param filename: name of file
    :type  filename: str
    
    :return: filename or path without suffix
    :rtype: str
    """
    try:
        if filename.find('.') != -1:
            filename = filename[: filename.rfind('.') ]     # remove ending
    except:
        pass  ## just in case there is no ending to start with...

    return filename


def stripFilename( filename ):
    """
    Return filename without path and without ending.

    :param filename: name of file
    :type  filename: str
    
    :return: base filename
    :rtype: str
    """
    name = osp.basename( filename )      # remove path
    try:
        if name.find('.') != -1:
            name = name[: name.rfind('.') ]     # remove ending
    except:
        pass  ## just in case there is no ending to start with...

    return name


def fileLength( filename ):
    """
    Count number of lines in a file.

    :param filename: name of file
    :type  filename: str
    
    :return: number of lines
    :rtype: int
    """
    p1 = subprocess.Popen( ['cat',filename], stdout=subprocess.PIPE )
    p2 = subprocess.Popen( ["wc", "-l"], stdin=p1.stdout,
                               stdout=subprocess.PIPE )
    return int(p2.communicate()[0])


def tempDir():
    """
    Get folder for temporary files - either from environment settings
    or '/tmp'

    :return: directort for temporary files
    :rtype: str    
    """
    if tempfile.tempdir is not None:
        return  tempfile.tempdir

    return osp.dirname( tempfile.mktemp() )


def file2dic( filename ):
    """
    Construct dictionary from file with key - value pairs (one per line).

    :param filename: name of file
    :type  filename: str
    
    :raise ToolsError: if file can't be parsed into dictionary
    :raise IOError: if file can't be opened
    """
    try:
        line = None
        result = {}
        for line in open( filename ):

            if '#' in line:
                line = line[ : line.index('#') ]
            line = line.strip()

            l = line.split()[1:]

            if len( l ) == 0 and len( line ) > 0:
                result[ line.split()[0] ] = ''
            if len( l ) == 1:
                result[ line.split()[0] ] = l[0]
            if len( l ) > 1:
                result[ line.split()[0] ] = l
    except:
        s = "Error parsing option file %s." % filename
        s += '\nLine: ' + str( line )
        s += '\n' + lastError()
        raise ToolsError( s )

    return result


def get_cmdDict(lst_cmd, dic_default):
    """
    Parse commandline options into dictionary of type C{ {<option> : <value>} }
    Options are recognised by a leading '-'.
    Error handling should be improved.
    
    Option C{ -x |file_name| } is interpreted as file with additional options.
    The key value pairs in lst_cmd replace key value pairs in the
    -x file and in dic_default.
    

    :param lst_cmd: list with the command line options::
                    e.g. ['-pdb', 'in1.pdb', 'in2.pdb', '-o', 'out.dat']
    :type  lst_cmd: [str]
    :param dic_default: dictionary with default options::
                        e.g. {'psf':'in.psf'}
    :type  dic_default: {str : str}

    :return: command dictionary::
             ala {'pdb':['in1.pdb', 'in2.pdb'], 'psf':'in.psf', 'o':'out.dat'}
    :rtype: {<option> : <value>}
    """
    dic_cmd = {}                     # create return dictionary
    try:

        for cmd in lst_cmd:
            if (cmd[0] == '-'):               # this entry is new option
                current_option = cmd[1:]      # take all but leading "-"
                dic_cmd[current_option] = ""  # make sure key exists even
                                              # w/o value
                counter = 0        # number of values for this option
            else:                  # this entry is value for latest option

                if counter < 1:
                    dic_cmd[current_option] = cmd

    # in case, several values follow after a "-xxx" option convert dictionary
    # entry into list and add all elements (until the next "-") to this list
                else:
                    if counter == 1:   # there is already a value assigned
    # convert to list
                        dic_cmd[current_option] = [dic_cmd[current_option]]
    # add value to list
                    dic_cmd[current_option] = dic_cmd[current_option] + [cmd]

                counter = counter + 1

    except (KeyError, UnboundLocalError) as why:
        errWriteln("Can't resolve command line options.\n \tError:"+str(why))

    ## get extra options from external file
    try:
        if 'x' in dic_cmd:
            d = file2dic( dic_cmd['x'] )
            d.update( dic_cmd )
            dic_cmd = d
    except IOError:
        errWriteln( "Error opening %s."% dic_cmd['x'] )
    except ToolsError as why:
        errWriteln( str(why) )

    ## fill in missing default values
    dic_default.update( dic_cmd )
    dic_cmd = dic_default

    return dic_cmd


def cmdDict( defaultDic={} ):
    """
    Convenience implementation of :class:`get_cmdDict`. Take command line options
    from sys.argv[1:] and convert them into dictionary.
    Example::
      '-o out.dat -in 1.pdb 2.pdb 3.pdb -d' will be converted to
      {'o':'out.dat', 'in': ['1.pdb', '2.pdb', '3.pdb'], 'd':'' }
      
    Option C{ -x |file_name| } is interpreted as file with additional options.
    
    :param defaultDic: dic with default values.
    :type  defaultDic: dic

    :return: command dictionary
    :rtype: dic
    """
    return get_cmdDict( sys.argv[1:], defaultDic )


def dump(this, filename, gzip = 0, mode = 'wb'):
    """
    Dump this::
      dump(this, filename, gzip = 0)
      Supports also '~' or '~user'.

    Note: Peter Schmidtke : gzip fixed, works now
    
    Written by Wolfgang Rieping.

    :param this: object to dump
    :type  this: any
    :param filename: name of file
    :type  filename: str
    :param gzip: gzip dumped object (default 0)
    :type  gzip: 1|0
    :param mode: file handle mode (default w)
    :type  mode: str
    """
    import biskit as B
    
    filename = osp.expanduser(filename)

    ## special case: do not slim PDBModels that are pickled to override
    ## their own source
    if isinstance( this, B.PDBModel ) and not this.forcePickle \
       and osp.samefile( str(this.source), filename ):
        this.saveAs( filename )

    else:

        if not mode in ['wb', 'w', 'a']:
            raise PickleError("mode has to be 'wb' (write binary) or 'a' (append)")

        if gzip:
            f = gzopen(filename, "wb")
        else:
            f = open(filename, mode)

        pickle.dump(this, f, 1)

        f.close()


def Dump( this, filename, gzip = 0, mode = 'w'):
    EHandler.warning('deprecated: tools.Dump has been renamed to tools.dump')
    return dump( this, filename, gzip=gzip, mode=mode )

def load(filename, gzip=0, encoding='ASCII'):
    """
    Load dumped object from file.

    Note: Peter Schmidtke : gzip fixed, works now
    
    Written by Wolfgang Rieping.

    :param filename: name of file
    :type  filename: str
    :param gzip: unzip dumped object (default 0)
    :type  gzip: 1|0
    :param encoding: optional encoding for pickle.load ['ASCII']
    :type encoding: str

    :return: loaded object
    :rtype: any

    :raise cPickle.UnpicklingError, if the pickle format is not recognized
    """
    filename = osp.expanduser(filename)

    try:
        if gzip :
            f=gzopen(filename,"rb")
        else :            
            f = open(filename, 'rb')

        objects = []
        eof = 0
        n = 0

        while not eof:
            try:
                this = pickle.load(f) #, encoding=encoding)
                objects.append(this)
                n += 1
            except EOFError:
                eof = 1

        f.close()

        if n == 1:
            return objects[0]
        else:
            return tuple(objects)

    except ValueError as why:
        raise PickleError('Python pickle %s is corrupted.' % filename)

def Load( filename, gzip=0 ):
    EHandler.warning('deprecated: tools.Load has been renamed to tools.load')
    return load( filename, gzip=gzip )



def packageRoot():
    """
    :return: absolute folder of the biskit python package.
    :rtype: str
    """
    import biskit
    return absfile( osp.split(biskit.__file__)[0] )

def projectRoot():
    """
    Root of biskit project. That's the folder **containing** the biskit python
    package (the parent folder of `packageRoot`).
    
    :return: absolute path of the root (parent) folder of current project::
             i.e. '/home/raik/py/biskitproject'
    :rtype: string
    """
    return absfile(osp.join( packageRoot(), '..' ))

def dataRoot():
    """
    Root of Biskit data directory (formerly 'biskit/Biskit/data').
    :return: absolute path
    :rtype: string    
    """
    return osp.join( packageRoot(), 'data' )

def testRoot( subfolder='' ):
    """
    Root of Biskit test directory.
    
    :param subfolder: str, optional sub-folder of test data folder

    :return: absolute path
    :rtype: string    
    """
    if subfolder and subfolder[0] == os.sep:
        subfolder = subfolder[1:]
    return os.path.join( packageRoot(), 'testdata', subfolder )


def isBinary( f ):
    """
    Check if file is a binary.
    
    :param f: path to existing file
    :type  f: str

    :return: condition
    :rtype: 1|0
    
    :raise OSError: if file doesn't exist
    """
    return os.access( f, os.X_OK )


def binExists( f ):
    """
    Check if binary with file name f exists.
    
    :param f: binary file name
    :type  f: str
    
    :return: True if binary file f is found in PATH and is executable
    :rtype: 1|0
    """
    if osp.exists( f ):
        return isBinary( f )

    for path in os.getenv( 'PATH' ).split(':') :

        full_path = osp.join( path, f )

        if osp.exists( full_path ) and isBinary( full_path ):
            return True

    return False


def absbinary( f ):
    """
    Absolute path of binary.
    
    :param f: binary file name
    :type  f: str
    
    :return: full path to existing binary
    :rtype: str
    
    :raise IOError: if an executable binary is not found in PATH
    """
    if osp.exists( f ) and isBinary( f ):
        return f

    for path in os.getenv( 'PATH' ).split(':') :

        full_path = osp.join( path, f )

        if osp.exists( full_path ) and isBinary( full_path ):
            return full_path

    raise IOError('binary %s not found.' % f)


def platformFolder( f ):
    """
    Get a platform-specific subfolder of f for platform-dependent imports.
    
    :param f: parent folder
    :type  f: str

    :return: path
    :rtype: str
    """
    import platform as P

    __version =  '.'.join( P.python_version().split('.')[:2] )

    r = 'py%s_%s' % (__version, P.machine() ) ## , P.architecture()[0] )

    r = os.path.join( f, r)

    return r


def sortString( s ):
    """
    Sort the letters of a string::
      sortString( str ) -> str with sorted letters
      
    :param s: string to be sorted
    :type  s: str
    
    :return: sorted string
    :rtype: str  
    """
    l = list(s)
    l.sort()
    return ''.join(l)


def string2Fname( s ):
    """
    Remove forbidden character from string so that it can be used as a
    filename.

    :param s: string
    :type  s: str

    :return: cleaned string
    :rtype: str     
    """
    forbidden = ['*', '?', '|', '/', ' ']
    replaceme = ['-', '-', '-', '-', '_']
    for i in range(0, len(forbidden)):
        s = s.replace( forbidden[i], replaceme[i] )
    return s


def toIntList( o ):
    """
    Convert single value or list of values into list of integers.
    
    :param o: value or list
    :type  o: int or [int]

    :return: list of integer
    :rtype: [int]
    """
    if type( o ) != type( [] ):
        o = [ o ]

    return list(map( int, o ))


def toIntArray( o ):
    """
    Convert single value or list of values to numpy array of int.
    
    :param o: value or list
    :type  o: int or [int]

    :return: array of integer
    :rtype: N0.array('i')    
    """
    if type( o ) == list or type( o ) == type( N0.array([])):
        return N0.array( map( int, o ) )

    return N0.array( [ int( o ) ] )


def toList( o ):
    """
    Make a list::
      toList(o) -> [o], or o,  if o is already a list
    
    :param o: value(s)
    :type  o: any or [any]

    :return: list
    :rtype: [any]
    """
    if type( o ) != type( [] ):
        return [ o ]
    return o


def toStr( o ):
    """
    Make a string from a list or interger.
    Stripping of any flanking witespaces.  
    
    :param o: value(s)
    :type  o: any or [any]

    :return: list
    :rtype: [any]
    """
    if type( o ) == type( 1 ):
        return str(o)
    
    if type( o ) == type( [] ):
        s = ''
        for item in o:
            s += str.strip( str(item) ) 
        return s
    
    return o


def toInt( o, default=None ):
    """
    Convert to intereg if possible::
      toInt(o) -> int, int(o) or default if o is impossible to convert.

    :param o: value
    :type  o: any
    :param default: value to return if conversion is impossible (default: None)
    :type  default: any
    
    :return: integer OR None 
    :rtype: int OR None     
    """
    if o is None or o == '':
        return default
    try:
        return int( o )
    except:
        return default


def hex2int( shex ):
    """
    Convert hex-code string into float number.
    :param s: hex-code, e.g. 'FF0B99'
    :type  s: str
    :return: float
    :rtype: float
    """
    shex = shex.replace('0x','')

    factors = [ 16**(i) for i in range(len(shex)) ]
    factors.reverse()
    factors = N0.array( factors )

    table   = dict( list(zip('0123456789abcdef',list(range(16)))) )

    components = [ table[s]*f for s,f in  zip( shex.lower(), factors ) ]

    return N0.sum( components )
    

def colorSpectrum( nColors, firstColor='FF0000', lastColor='FF00FF' ):
    """
    Creates a list of 'nColors' colors for biggles starting at
    'firstColor' ending at 'lastColor'
    Examples::
     free spectrum red FF0000 to green 00FF00
     bound spectrum cyan 00FFFF to magenta FF00FF

    :param nColors: number of colors to create
    :type  nColors: int
    :param firstColor: first color in hex format (default: FF0000)
    :type  firstColor: str
    :param lastColor: last color in hex format (default: FF00FF)
    :type  lastColor: str
    
    :return: list of colors
    :rtype: [int]
    """
    spec = []
    cmd = dataRoot() + '/spectrum.pl ' +str(nColors) +\
                    ' ' + str(firstColor) + ' ' +  str(lastColor)
    out = os.popen( cmd ).readlines()
    
    if not out:
        raise IOError('color generation command failed: %r' % cmd)

    for s in out:    
        spec += [ hex2int( str( str.strip( s  ) ) ) ]

    return spec


def hexColors( nColors, firstColor='FF0000', lastColor='FF00FF' ):
    """
    Creates a list of 'nColors' colors for PyMol starting at
    'firstColor' ending at 'lastColor'
    Examples::
     free spectrum red FF0000 to green 00FF00
     bound spectrum cyan 00FFFF to magenta FF00FF
     
    :param nColors: number of colors to create
    :type  nColors: int
    :param firstColor: first color in hex format (default: FF0000)
    :type  firstColor: str
    :param lastColor: last color in hex format (default: FF00FF)
    :type  lastColor: str

    :return: list of  hex colors
    :rtype: [ str ]
    """
    spec = []
    cmd = dataRoot() + '/spectrum.pl ' +str(nColors) +\
                    ' ' + str(firstColor) + ' ' +  str(lastColor)
    out = os.popen( cmd ).readlines()

    if not out:
        raise IOError('color generation command failed: %r' % cmd)

    for s in out:    
        spec += [ '0x' + str( str.strip( s  ) ) ]

    return spec


def rgb2hex( rgbColor ):
    """
    convert rgb color into 8 bit hex rgb color::
      [ 1.0, 0.0, 1.0, ] -> 'FF00FF'

    :param rgbColor: RGB-color e.g. [ 1.0, 0.0, 1.0, ]
    :type rgbColor : [float]
    
    :return: hex colors
    :rtype: str 
    """
    hexRgb = ''
    for i in range(0,3):        
        component = hex( int( rgbColor[i]*255 ) )[2:]

        if len(component) == 1:
            hexRgb += '0' + component
        else:
            hexRgb +=  component

    return hexRgb


def hex2rgb( hexColor, str=0 ):
    """
    convert 8 bit hex rgb color into  rgb color  ::
       'FF00FF' -> [ 1.0, 0.0, 1.0, ] 

    :param hexColor: HEX-color e.g. 'FF00FF'
    :type  hexColor: str
    :param str: return rgb colors as a tring (i.e for PyMol)
    :type  str: 1|0
    
    :return: rgb colors
    :rtype: [float]
    """
    rgb = []
    if hexColor[:2] == '0x':
        hexColor = hexColor[2:]
    
    for i in range(0,6,2):
        rgb += [ int(hexColor[i:i+2], 16)/255.0 ]

    if str:
        rgb_str= '[ %.3f, %.3f, %.3f ]'%(rgb[0], rgb[1], rgb[2])
        
        return rgb_str
    
    return rgb


def dateString():
    """
    :return: DD/MM/YYYY
    :rtype: str
    """
    t = localtime()
    return '%02i/%02i/%i' % (t[2],t[1],t[0] )


def dateSortString():
    """
    :return: YYYY/MM/DD:hh:mm.ss.ms
    :rtype: 
    """
    t = localtime()
    return "%i/%02i/%02i:%02i.%02i.%02i" % (t[0],t[1],t[2],t[3],t[4],t[5])


def tryRemove(f, verbose=0, tree=0, wildcard=0 ):
    """
    Remove file or folder::
     remove(f [,verbose=0, tree=0]), remove if possible, otherwise do nothing

    :param f: file path
    :type  f: str
    :param verbose: report failure (default 0)
    :type  verbose: 0|1
    :param tree: remove whole folder (default 0)
    :type  tree: 0|1
    :param wildcard: filename contains wildcards (default 0)
    :type  wildcard: 0|1

    :return: 1 if file was removed
    :rtype: 1|0
    """
    try:
        if osp.isdir(f):
            if tree:
                shutil.rmtree( f, ignore_errors=1 )
            else:
                errWriteln('%r is directory - not removed.' % f)
        else:
            if wildcard:
                l = glob.glob( f )
                for i in l:
                    os.remove( i )
            else:
                os.remove( f )
        return 1
    except:
        if verbose: errWriteln( 'Warning: Cannot remove %r.' % f )
        return 0

def backup( fname, suffix='~' ):
    """
    Create backup of file if it already exists.
    :param fname: file name
    :type  fname: str
    :param suffix: suffix to add to backup file name ['~']
    :type  suffix: str

    :return: True if backup was created, False otherwise
    :rtype: bool
    """
    fname = absfile( fname )
    
    if os.path.exists( fname ):
        os.rename( fname, fname + suffix )
        return True
    return False


def ensure( v, t, allowed=[], forbidden=[] ):
    """
    Check type of a variable
    
    :param v: variable to test
    :type  v: variable
    :param t: required type
    :type  t: str
    :param allowed: list of additional values allowed for v {default: []}
    :type  allowed: [str]
    
    :raise TypeError: if invalid
    """
    if allowed:
        allowed = toList( allowed )
        if len( allowed ) > 0 and v in allowed:
            return

    if not isinstance(v, t):
        raise TypeError('looked for %s but found %s' % (str(t),str(v)[:20]))

    if forbidden and v in forbidden:
        raise TypeError('value %s is not allowed.' % (str(v)[:20]))


def clipStr( s, length, suffix='..', expandtabs=1 ):
    """
    Shorten string from end and replace the last characters with suffix::
      clipStr( str, length ) -> str, with len( str ) <= length

    :param s: original string
    :type  s: str
    :param length: desired length
    :type  length: int
    :param suffix: suffix (default: ..)
    :type  suffix: str
    
    :return: shortend string
    :rtype: str
    """
    if expandtabs:
        s = s.expandtabs()
        
    if len(s) > length:
        s = s[:(length - len(suffix))] + suffix
    return s


def info( item, short=1 ):
    """
    ::
      info( item, short=1) -> Print useful information about item.

    :param item: query item
    :type  item: item
    :param short: short version (default: 1)
    :type  short: 1|0
    """
    ## quick and dirty ##
    if hasattr(item, '__name__'):
        print("NAME:    ", item.__name__)
    if hasattr(item, '__class__'):
        print("CLASS:   ", item.__class__.__name__)
    print("ID:      ", id(item))
    print("TYPE:    ", type(item))
    print("VALUE:   ", repr(item))
    print("CALLABLE:", end=' ')
    if isinstance(item, collections.abc.Callable):
        print("Yes")
    else:
        print("No")
    if hasattr(item, '__doc__'):
        doc = getattr(item, '__doc__')
        if doc:
            doc = doc.strip()   # Remove leading/trailing whitespace.
            if short:
                doc = doc.split('\n')[0]
            print("DOC:    ", '\n\t' * (not short), doc)

    print("\nMETHODS")
    methods = [ getattr( item, m ) for m in dir( item )
                if isinstance( getattr( item, m ), collections.abc.Callable) ]

    for m in methods:
        doc = getattr(m, '__doc__', '')
        if doc:
            doc = str(doc).strip()
            if short:
                doc = str(doc).split('\n')[0]
        else:
            doc = ''
        s = "%-15s: " % (getattr(m,'__name__','?')) + '\n\t'*(not short) + doc
        if short:
            s = clipStr( s, 79 )
        print(s)

    if hasattr( item, '__dict__'):

        print("\nFIELDS")
        for k, v in list(item.__dict__.items()):
            s = "%-15s: %s" % (k, str(v).strip() )
            print(clipStr( s, 79 ))


class PseudoClass(object):
    """
    Empty class that raises an ImportError upon creation.
    """
    def __new__(cls, *args, **kw):
        raise ImportError('Class %r is not available because of missing modules: %r' \
              % (cls.__name__, str(cls.error)))
    
def tryImport( module, cls, old_as=None, namespace=None ):
    """
    Try to import a class from a module. If that fails, 'import' a
    default class of the same name that raises an exception when used.
    
    :param module: name of the module
    :type  module: str
    :param cls  : name of the class
    :type  cls  : str
    :param namespace: namespace for the import [default: globals() ]
    :type  namespace: dict
    
    :return: True if import succeeded, False otherwise
    :rtype: bool
    """
    old_as = old_as or cls
    g = namespace or globals()
    try:
        exec('from %s import %s as %s' % (module, cls, old_as), g)
        return True
    
    except ImportError as e:

        Cls = type( cls,(PseudoClass,),{'error':e} )
        g.update( {old_as: Cls} )

    return False

def tryImportModule( module, old_as=None, namespace=None ):
    """
    Try to import a module. If that fails, 'import' a dummy module
    of the same name.

    NOTE: as of python 3, the namespace returned by globals() cannot
    reach out of the tools module -- this method thus always needs to be
    called like this:

    >>> tryImportModule( 'numpy', namespace=globals() )
    
    :param module: name of the module
    :type  module: str
    :param namespace: namespace for the import [default: globals() ]
    :type  namespace: dict
    
    :return: True if import succeeded, False otherwise
    :rtype: bool
    """
    old_as = old_as or module
    g = namespace or globals()
    try:
        exec('import %s as %s' % (module, old_as), g)
        return True
    
    except ImportError as e:

        m = types.ModuleType( old_as, doc='Pseudo module. Import of real one failed.' )
        m.error = str(e)
        
        g.update( {old_as: m} )

    return False

def gzopen( fname, mode='rt' ):
    """
    Open a normal or a gzipped file.
    :param fname: file name (can contain ~, .. etc.)
    :type  fname: str
    :param mode: read/write mode ['r']
    :type  mode: str
    
    :return: file handle
    :rtype: file or GZipFile
    """
    fname = absfile( fname )

    if fname[-2:] == 'gz':
        return gzip.open( fname, mode )

    return open( fname, mode )


def profile( s ):
    """
    Profile the given code fragment and report time-consuming method calls.
    :param s: python code fragment, example: 'm = PDBModel("3tgi")'
    :type  s: str
    """    
    try:
        import cProfile as profile
        fout = tempDir() + '/profiling.out'

        profile.run( s, fout )

        ## Analyzing
        import pstats
        p = pstats.Stats(fout)

        ## long steps and methods calling them
        p.sort_stats('cumulative').print_stats(20)
        p.print_callers(0.0)

        tryRemove( fout )

    except ImportError as why:
        raise ToolsError('Python profiling modules are not installed.')


#############
##  TESTING        
#############
import biskit.test as BT
       
class Test(BT.BiskitTest):
    """Test case"""

    def prepare(self):
        import tempfile
        self.f_pickle = tempfile.mktemp('.pickle')

    def cleanUp( self ):
        tryRemove( self.f_pickle )

    def test_error_reporting( self ):
        try:
            i = 1/0
        except:
            if self.local:
                print('\nTest error trace:\n' + lastErrorTrace())
            self.assertTrue( lastErrorTrace() != '' )

    def test_ensure_types(self):
        self.assertIsNone( ensure( 'teststr', str ) )
        
        with self.assertRaises(TypeError):
            ensure( 'teststr', int )
    
    def test_filehandling(self):
        f = absfile('~')
        self.assertTrue(osp.exists(f))

        self.fpdb = testRoot( '/rec/1A2P.pdb' )
        if self.VERBOSITY > 2:
            self.log.write('path: ' + self.fpdb)
        self.assertTrue( osp.exists(self.fpdb) )

    def test_isBinary(self):
        self.assertEqual( isBinary('/usr/bin/nice'), 1)
    
    EXPECT_colors = [16711680, 15493120, 11259136, 5570304, 458570, 63925, 
                     44031, 3298559, 9901055, 16711935]
    EXPECT_hexcolors = ['0xFF0000', '0xEC6800', '0xABCD00', '0x54FF00', 
                        '0x06FF4A', '0x00F9B5', '0x00ABFF', '0x3254FF', 
                        '0x9713FF', '0xFF00FF']
    def test_colorSpectrum(self):
        colors = colorSpectrum(10)
        self.assertAlmostEqual(colors, self.EXPECT_colors)
        
        hexcolors = hexColors(10)
        self.assertEqual(hexcolors, self.EXPECT_hexcolors)
    
    def test_pickling(self):
        d1 = {'a':10, 'b':10.0, 'c':'test'}
        dump(d1, self.f_pickle )
        d2 = load(self.f_pickle)
        
        self.assertEqual(d1, d2)

if __name__ == '__main__':
    
    BT.localTest()

