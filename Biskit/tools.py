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
## last $Author$
## last $Date$
## $Revision$

"""
simply -  everyday tools
"""

import sys, string
import os.path as osp
import shutil
import operator
import os, cPickle  ## for Load and Dump
import tempfile
import traceback
from inspect import getframeinfo
from time import localtime
import Numeric
import types
import glob
import subprocess

class ToolsError( Exception ):
    pass


def errWriteln(s):
    """
    print s to standard error with line feed.
    
    @param s: string
    @type  s: str
    """
    sys.stderr.write(s+'\n')
    sys.stderr.flush()


def errWrite(s):
    """
    print s to standard error.

    @param s: string
    @type  s: str    
    """
    sys.stderr.write(s)
    sys.stderr.flush()


def flushPrint(s):
    """
    print s without line break and flush standard out.

    @param s: string
    @type  s: str    
    """
    sys.stdout.write(s)
    sys.stdout.flush()


def lastError():
    """
    Collect type and line of last exception.
    
    @return: '<ExceptionType> in line <lineNumber>:<Exception arguments>'
    @rtype: String
    """
    try:
        trace = sys.exc_info()[2]
        why = sys.exc_info()[1]
        try:
            why = sys.exc_info()[1].args
        except:
            pass
        file = getframeinfo( trace.tb_frame )[0]

        result = "%s in %s line %i:\n\t%s." % ( str(sys.exc_type),
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


def dictAdd( dic, key, value ):
    """
    Add value to dic, create list, if dic has already value in key.

    @param key: dictionary key
    @type  key: str
    @param value: value
    @type  value: any
    """
    if key in dic:
        old = dic[key]

        if type( old ) != list and value != old:
            dic[ key ] = [ old ] + [ value ]
        else:
            if type( old ) == list and value not in old:
                dic[ key ] = old + [ value ]

    else:
        dic[key] = value


def absfile( filename, resolveLinks=1 ):
    """
    Get absolute file path::
      - expand ~ to user home, change
      - expand ../../ to absolute path
      - resolve links

    @param filename: name of file
    @type  filename: str
    @param resolveLinks: eliminate any symbolic links (default: 1)
    @type  resolveLinks: 1|0
    
    @return: absolute path or filename
    @rtype: string
    """
    if not filename:
        return filename
    r = osp.abspath( osp.expanduser( filename ) )
    if resolveLinks:
        r = osp.realpath( r )
    return r


def homefile( filename, otherUser=1, ownCopy=1 ):
    """
    Relativize a file name to ~ or, if it is in another user's home,
    to ~otheruser or, if it is in nobody's home, to / .
    
    L{splithome()} is used to also guess home directories of other users.
    
    @param filename: name of file
    @type  filename: str
    @param otherUser: look also in other user's home directories (default 1)
    @type  otherUser: 1|0
    @param ownCopy: replace alien path by path into own home directory if
                    possible, e.g. ~other/data/x is replaced
                    by ~/data/x if there is such a file. (default 1) Careful!
    @type  ownCopy: 1|0
                
    @return: path or filename
    @rtype: str
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

    @param filename: name of file
    @type  filename: str
    
    @return: home folder of some user, remaining path relative to home
    @rtype: (str, str)
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


def stripSuffix( filename ):
    """
    Return file name without ending.

    @param filename: name of file
    @type  filename: str
    
    @return: filename or path without suffix
    @rtype: str
    """
    try:
        if filename.find('.') <> -1:
            filename = filename[: filename.rfind('.') ]     # remove ending
    except:
        pass  ## just in case there is no ending to start with...

    return filename


def stripFilename( filename ):
    """
    Return filename without path and without ending.

    @param filename: name of file
    @type  filename: str
    
    @return: base filename
    @rtype: str
    """
    name = osp.basename( filename )      # remove path
    try:
        if name.find('.') <> -1:
            name = name[: name.rfind('.') ]     # remove ending
    except:
        pass  ## just in case there is no ending to start with...

    return name


def fileLength( filename ):
    """
    Count number of lines in a file.

    @param filename: name of file
    @type  filename: str
    
    @return: number of lines
    @rtype: int
    """
    p1 = subprocess.Popen( ['cat',filename], stdout=subprocess.PIPE )
    p2 = subprocess.Popen( ["wc", "-l"], stdin=p1.stdout,
                               stdout=subprocess.PIPE )
    return int(p2.communicate()[0])


def tempDir():
    """
    Get folder for temporary files - either from environment settings
    or '/tmp'

    @return: directort for temporary files
    @rtype: str    
    """
    if tempfile.tempdir != None:
        return  tempfile.tempdir

    return osp.dirname( tempfile.mktemp() )


def file2dic( filename ):
    """
    Construct dictionary from file with key - value pairs (one per line).

    @param filename: name of file
    @type  filename: str
    
    @raise ToolsError: if file can't be parsed into dictionary
    @raise IOError: if file can't be opened
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
        s = "Error parsing option file %s." % fname
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
    

    @param lst_cmd: list with the command line options::
                    e.g. ['-pdb', 'in1.pdb', 'in2.pdb', '-o', 'out.dat']
    @type  lst_cmd: [str]
    @param dic_default: dictionary with default options::
                        e.g. {'psf':'in.psf'}
    @type  dic_default: {str : str}

    @return: command dictionary::
             ala {'pdb':['in1.pdb', 'in2.pdb'], 'psf':'in.psf', 'o':'out.dat'}
    @rtype: {<option> : <value>}
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

    except (KeyError, UnboundLocalError), why:
        errWriteln("Can't resolve command line options.\n \tError:"+str(why))

    ## get extra options from external file
    try:
        if dic_cmd.has_key('x'):
            d = file2dic( dic_cmd['x'] )
            d.update( dic_cmd )
            dic_cmd = d
    except IOError:
        errWriteln( "Error opening %s."% dic_cmd['x'] )
    except ToolsError, why:
        errWriteln( str(why) )

    ## fill in missing default values
    dic_default.update( dic_cmd )
    dic_cmd = dic_default

    return dic_cmd


def cmdDict( defaultDic={} ):
    """
    Convenience implementation of L{get_cmdDict}. Take command line options
    from sys.argv[1:] and convert them into dictionary.
    Example::
      '-o out.dat -in 1.pdb 2.pdb 3.pdb -d' will be converted to
      {'o':'out.dat', 'in': ['1.pdb', '2.pdb', '3.pdb'], 'd':'' }
      
    Option C{ -x |file_name| } is interpreted as file with additional options.
    
    @param defaultDic: dic with default values.
    @type  defaultDic: dic

    @return: command dictionary
    @rtype: dic
    """
    return get_cmdDict( sys.argv[1:], defaultDic )


def Dump(this, filename, gzip = 0, mode = 'w'):
    """
    Dump this::
      Dump(this, filename, gzip = 0)
      Supports also '~' or '~user'.

    @author: Wolfgang Rieping
    @note gzip currently doesn't work.

    @param this: object to dump
    @type  this: any
    @param filename: name of file
    @type  filename: str
    @param gzip: gzip dumped object (default 0)
    @type  gzip: 1|0
    @param mode: file handle mode (default w)
    @type  mode: str
    """
    filename = osp.expanduser(filename)

    if not mode in ['w', 'a']:
        raise "mode has to be 'w' (write) or 'a' (append)"

##     if gzip:
##         import gzip
##         f = gzip.GzipFile(filename, mode)
##     else:
    f = open(filename, mode)

    cPickle.dump(this, f, 1)

    f.close()


def Load(filename, gzip = 0):
    """
    Load dumped object from file.

    @author: Wolfgang Rieping

    @param filename: name of file
    @type  filename: str
    @param gzip: unzip dumped object (default 0)
    @type  gzip: 1|0

    @return: loaded object
    @rtype: any
    """
    filename = osp.expanduser(filename)

    f = open(filename)

    objects = []

    eof = 0
    n = 0

    while not eof:
        try:
            this = cPickle.load(f)
            objects.append(this)
            n += 1
        except EOFError:
            eof = 1

    f.close()

    if n == 1:
        return objects[0]
    else:
        return tuple(objects)


## obsolete
def getOnDemand( attr, dumpIt=1):
    """
    Return attr either unpickled, freshly calculated or unchanged.
    If attr is no tuple or anything else goes wrong it is returned unchanged
    
    @param attr: tuple of string and function::
                 attr[0] - name of existing or non-existing file or ''
                 attr[1] - function to get result if it can't be unpickled
    @type  attr: (str, function)
    @param dumpIt: try pickling return value to attr[0](if valid file)
                   (default 1)
    @type  dumpIt: 0|1

    @return: attr (unchanged or object unpickeled from file)
    @rtype: (str, function)
    """
    try:
        if type( attr ) == type( ('','') ) and type( attr[0] ) == type( '' ):

            fname = attr[0]
            function = attr[1]

            ## file exists, try unpickling from it
            if type(fname)==type('') and osp.exists( attr[0] ):
                return Load( fname )

            ## pretend function exists, try to calculate return value
            result = function()

            ## pickle result if file would end up in a valid path
            if osp.exists( osp.dirname( fname )) and dumpIt:
                Dump( result, fname )

            return result

    except:
        print lastError()

    ## return attr directly if above failed
    return attr


def projectRoot():
    """
    Root of biskit project.
    
    @return: absolute path of the root of current project::
             i.e. '/home/Bis/raik/biskit'
    @rtype: string
    """
    ## import this module
    import tools
    ## get location of this module
    f = absfile(tools.__file__)
    ## extract path and assume it is 'project_root/Biskit'
    f = osp.split( f )[0] + '/../'
    return absfile( f )


def testRoot():
    """
    Root of Biskit test directory.

    @return: absolute path
    @rtype: string    
    """
    return projectRoot() + '/test'


def isBinary( f ):
    """
    Check if file is a binary.
    
    @param f: path to existing file
    @type  f: str

    @return: condition
    @rtype: 1|0
    
    @raise OSError: if file doesn't exist
    """
    st = os.stat( f )
    mode = st[0]                         ## permissions
    mode = int( oct( mode & 0777 )[1] )  ## user permissions as octal value

    return ( operator.and_( mode, 1 ) == 1 )


def binExists( f ):
    """
    Check if binary with file name f exists.
    
    @param f: binary file name
    @type  f: str
    
    @return: True if binary file f is found in PATH and is executable
    @rtype: 1|0
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
    
    @param f: binary file name
    @type  f: str
    
    @return: full path to existing binary
    @rtype: str
    
    @raise IOError: if an executable binary is not found in PATH
    """
    if osp.exists( f ) and isBinary( f ):
        return f

    for path in os.getenv( 'PATH' ).split(':') :

        full_path = osp.join( path, f )

        if osp.exists( full_path ) and isBinary( full_path ):
            return full_path

    raise IOError, 'binary %s not found.' % f


def platformFolder( f ):
    """
    Get a platform-specific subfolder of f for platform-dependent imports.
    
    @param f: parent folder
    @type  f: str

    @return: path
    @rtype: str
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
      
    @param s: string to be sorted
    @type  s: str
    
    @return: sorted string
    @rtype: str  
    """
    l = list(s)
    l.sort()
    return ''.join(l)


def string2Fname( s ):
    """
    Remove forbidden character from string so that it can be used as a
    filename.

    @param s: string
    @type  s: str

    @return: cleaned string
    @rtype: str     
    """
    forbidden = ['*', '?', '|', '/', ' ']
    replaceme = ['-', '-', '-', '-', '_']
    for i in range(0, len(forbidden)):
        s = s.replace( forbidden[i], replaceme[i] )
    return s


def toIntList( o ):
    """
    Convert single value or list of values into list of integers.
    
    @param o: value or list
    @type  o: int or [int]

    @return: list of integer
    @rtype: [int]
    """
    if type( o ) != type( [] ):
        o = [ o ]

    return map( int, o )


def toIntArray( o ):
    """
    Convert single value or list of values to Numeric array of int.
    
    @param o: value or list
    @type  o: int or [int]

    @return: array of integer
    @rtype: Numeric.array('i')    
    """
    if type( o ) == list or type( o ) == type( Numeric.array([])):
        return Numeric.array( map( int, o ) )

    return Numeric.array( [ int( o ) ] )


def toList( o ):
    """
    Make a list::
      toList(o) -> [o], or o,  if o is already a list
    
    @param o: value(s)
    @type  o: any or [any]

    @return: list
    @rtype: [any]
    """
    if type( o ) != type( [] ):
        return [ o ]
    return o


def toInt( o, default=None ):
    """
    Convert to intereg if possible::
      toInt(o) -> int, int(o) or default if o is impossible to convert.

    @param o: value
    @type  o: any
    @param default: value to return if conversion is impossible (default: None)
    @type  default: any
    
    @return: integer OR None 
    @rtype: int OR None     
    """
    if o == None or o == '':
        return default
    try:
        return int( o )
    except:
        return default


def colorSpectrum( nColors, firstColor='FF0000', lastColor='FF00FF' ):
    """
    Creates a list of 'nColors' colors for biggles starting at
    'firstColor' ending at 'lastColor'
    Examples::
     free spectrum red FF0000 to green 00FF00
     bound spectrum cyan 00FFFF to magenta FF00FF

    @param nColors: number of colors to create
    @type  nColors: int
    @param firstColor: first color in hex format (default: FF0000)
    @type  firstColor: str
    @param lastColor: last color in hex format (default: FF00FF)
    @type  lastColor: str
    
    @return: list of colors
    @rtype: [int]
    """
    spec = []
    out = os.popen( projectRoot() + '/external/spectrum.pl ' +str(nColors) +
                    ' ' + str(firstColor) + ' ' +  str(lastColor) ).readlines()

    for s in out:    
        spec += [ int( float( '0x' + str( string.strip( s  ) ) ) ) ]

    return spec


def hexColors( nColors, firstColor='FF0000', lastColor='FF00FF' ):
    """
    Creates a list of 'nColors' colors for PyMol starting at
    'firstColor' ending at 'lastColor'
    Examples::
     free spectrum red FF0000 to green 00FF00
     bound spectrum cyan 00FFFF to magenta FF00FF
     
    @param nColors: number of colors to create
    @type  nColors: int
    @param firstColor: first color in hex format (default: FF0000)
    @type  firstColor: str
    @param lastColor: last color in hex format (default: FF00FF)
    @type  lastColor: str

    @return: list of  hex colors
    @rtype: [ str ]
    """
    spec = []
    out = os.popen( projectRoot() + '/external/spectrum.pl ' +str(nColors) +
                    ' ' + str(firstColor) + ' ' +  str(lastColor) ).readlines()

    for s in out:    
        spec += [ '0x' + str( string.strip( s  ) ) ]

    return spec


def rgb2hex( rgbColor ):
    """
    convert rgb color into 8 bit hex rgb color::
      [ 1.0, 0.0, 1.0, ] -> 'FF00FF'

    @param rgbColor: RGB-color e.g. [ 1.0, 0.0, 1.0, ]
    @type rgbColor : [float]
    
    @return: hex colors
    @rtype: str 
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

    @param hexColor: HEX-color e.g. 'FF00FF'
    @type  hexColor: str
    @param str: return rgb colors as a tring (i.e for PyMol)
    @type  str: 1|0
    
    @return: rgb colors
    @rtype: [float]
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
    @return: DD/MM/YYYY
    @rtype: str
    """
    t = localtime()
    return '%02i/%02i/%i' % (t[2],t[1],t[0] )


def dateSortString():
    """
    @return: YYYY/MM/DD:hh:mm.ss.ms
    @rtype: 
    """
    t = localtime()
    return "%i/%02i/%02i:%02i.%02i.%02i" % (t[0],t[1],t[2],t[3],t[4],t[5])


def tryRemove(f, verbose=0, tree=0, wildcard=0 ):
    """
    Remove file or folder::
     remove(f [,verbose=0, tree=0]), remove if possible, otherwise do nothing

    @param f: file path
    @type  f: str
    @param verbose: report failure (default 0)
    @type  verbose: 0|1
    @param tree: remove whole folder (default 0)
    @type  tree: 0|1
    @param wildcard: filename contains wildcards (default 0)
    @type  wildcard: 0|1

    @return: 1 if file was removed
    @rtype: 1|0
    """
    try:
        if osp.isdir(f):
            if tree:
                shutil.rmtree( f, ignore_errors=1 )
            else:
                errWriteln('%s is directory - not removed.')
        else:
            if wildcard:
                l = glob.glob( f )
                for i in l:
                    os.remove( i )
            else:
                os.remove( f )
        return 1
    except:
        if verbose: errWriteln( 'Warning: Cannot remove %s.' % str(f) )
        return 0


def ensure( v, t, allowed=[], forbidden=[] ):
    """
    Check type of a variable
    
    @param v: variable to test
    @type  v: variable
    @param t: required type
    @type  t: str
    @param allowed: list of additional values allowed for v {default: []}
    @type  allowed: [str]
    
    @raise TypeError: if invalid
    """
    if allowed:
        allowed = toList( allowed )
        if len( allowed ) > 0 and v in allowed:
            return

    if not isinstance(v, t):
        raise TypeError, 'looked for %s but found %s' % (str(t),str(v)[:20])

    if forbidden and v in forbidden:
        raise TypeError, 'value %s is not allowed.' % (str(v)[:20])


def clipStr( s, length, suffix='..' ):
    """
    Shorten string from end and replace the last characters with suffix::
      clipStr( str, length ) -> str, with len( str ) <= length

    @param s: original string
    @type  s: str
    @param length: desired length
    @type  length: int
    @param suffix: suffix (default: ..)
    @type  suffix: str
    
    @return: shortend string
    @rtype: str
    """
    if len( s ) > length:
        s = s[:(length - len(suffix))] + suffix
    return s


def info( item, short=1 ):
    """
    Print info about ithem::
      info( item, short=1) -> Print useful information about item.

    @param item: query item
    @type  item: item
    @param short: short version (default: 1)
    @type  short: 1|0
    """
    ## quick and dirty ##
    if hasattr(item, '__name__'):
        print "NAME:    ", item.__name__
    if hasattr(item, '__class__'):
        print "CLASS:   ", item.__class__.__name__
    print "ID:      ", id(item)
    print "TYPE:    ", type(item)
    print "VALUE:   ", repr(item)
    print "CALLABLE:",
    if callable(item):
        print "Yes"
    else:
        print "No"
    if hasattr(item, '__doc__'):
        doc = getattr(item, '__doc__')
        if doc:
            doc = doc.strip()   # Remove leading/trailing whitespace.
            if short:
                doc = doc.split('\n')[0]
            print "DOC:    ", '\n\t' * (not short), doc

    print "\nMETHODS"
    methods = [ getattr( item, m ) for m in dir( item )
                if callable( getattr( item, m ) ) ]

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
        print s

    if hasattr( item, '__dict__'):

        print "\nFIELDS"
        for k, v in item.__dict__.items():
            s = "%-15s: %s" % (k, str(v).strip() )
            print clipStr( s, 79 )


###################################################
# main
###################################################
if __name__ == '__main__':

    print "TEST Exception:"
    try:
        i = 1/0
    except:
        print lastErrorTrace()
