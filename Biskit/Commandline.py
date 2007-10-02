## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2007 Raik Gruenberg & Johan Leckner
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
## last $Author: graik $
## last $Date: 2007-09-30 19:58:21 +0200 (Sun, 30 Sep 2007) $
## $Revision: 536 $
"""
parse and organize commandline options
"""

import sys
import Biskit.tools as T
import Biskit as B

class CommandlineError( B.BiskitError ):
    pass

class Commandline(dict):
    """
    Parse and organize commandline options.
    
    """
    
    def __init__( self, **kw ):
        self.defaults = {}
        self.types = {}
        self.docs = {}
        self.test = {}
        
        self.__setDefaults( **kw )

    def __parseDefault( self, k, v ):

        tdefault = type(v) if v is not None else str #: default type
        vdefault = v       #: default default value

        if type(v) is type:
            tdefault = v
            vdefault = None
            if tdefault is bool:
                vdefault = False

        elif type(v) is tuple:
            assert type(v[0]) is type, '%r is not a type' % v[0]
            assert len(v) == 2, 'tuple must be (type, defaultvalue)'
            tdefault = v[0]
            vdefault = v[1]
        
        elif type(v) is list:
            r = [ self.__parseDefault( k, x ) for x in v ]
            vdefault = [ x[1] for x in r ]
            tdefault = [ r[0][0] ]
        
        return tdefault, vdefault
    
    def __setDefaults( self, **kw ):
        """Called automatically by __init__"""

        for k, v in kw.items():
            self.types[k], self.defaults[k] = self.__parseDefault( k, v )

    
    def setTest( self, **testparam ):
        """Provide arguments for a unittest."""
        self.test = testparam

    def setDocs( self, **docs ):
        """
        Provide a short description of command line arguments.
        @param docs - key=value pairs (argument='description')
        @type docs  - { str:str }
        """
        self.docs = docs
    
    def weakupdate( self, d ):
        """
        Copy over values from d but don't override existing values
        @param d: dictionary with additional or default items
        @type  d: dict
        """
        for k, v in d.items():
            if not k in self:
                self[k] = v
    
    def isvalid( self, k, v ):
        """
        @param k: argument key
        @type  k: str
        @param v: value
        @type  v: any
        @return: True if v is a valid value for argument k
        """
        t = self.types.get( k, str )

        if type( t ) is list:
            t = t[0]
        assert type( t ) is type

        if type( v ) is list:
            r = True        
            for i in v:
                r = r and self.isvalid( k, i )
            return r

        return type( v ) is self.types.get( k, str )
    
    def validate( self ):
        """
        Validate all items of this dictionary (call .isvalid on all of them).
        
        Note: If an argument has a list of values, each item is validated
              separately.
        
        @raise CommandLineError, on the first invalid item
        """
        for k, v in self.items():
            if not self.isvalid( k, v ):
                raise CommandLineError('value %r is not valid for argument %r'\
                                       % (k,v))

    def __typecast( self ):
        for k, v in self.items():
            t = self.types.get( k, type(v) )
            try:
                
                if type(t) is list:
                    t = t[0]
                    if type(v) is not list:
                        v = [ v ]
                    self[k] = [ t(i) for i in v ]

                else:
                    self[ k ] = t( v )
            except ValueError, error:
                raise CommandLineError("argument %r (%r) doesn't fit %r"\
                                       % (k, v, t) )
            except TypeError, error:
                raise CommandLineError("argument %r (%r) doesn't fit %r"\
                                       % (k, v, t) )

    def __file2dic( self, f ):
        """
        @param f: input file, table of 'argument /t value'
        @type f: str
        @raise T.ToolsError, if anything goes wrong
        """
        return T.file2dic( f )
    
    def updateFromFile( self, f ):
        """
        Parse 'commandline arguments' from file.
        """
        ## get extra options from external file but don't override anything
        try:
            d = self.__file2dic( f )
            self.weakupdate( d )
        except IOError:
            raise CommandLineError( "Error opening %r."% f )
        except T.ToolsError, why:
            raise CommandLineError( str(why) )

    def __parse( self, argv ):
        """
        Parse commandline options into dictionary of type 
        C{ {<option> : <value>} }
        Options are recognised by a leading '-'.
        
        Option C{ -x |file_name| } is interpreted as file with additional options.
        The key value pairs in lst_cmd replace key value pairs in the
        -x file and in dic_default.
    
        @param lst_cmd: list with the command line options::
                        e.g. ['-pdb', 'in1.pdb', 'in2.pdb', '-o', 'out.dat']
        @type  lst_cmd: [str]
    
        @return: command dictionary::
                 ala {'pdb':['in1.pdb', 'in2.pdb'], 'psf':'in.psf', 'o':'out.dat'}
        @rtype: {<option> : <value>}
        """
        argv = argv[1:]
        try:
    
            for cmd in argv:
                if (cmd[0] == '-'):            # this entry is new option
                    current_option = cmd[1:]    # take all but leading "-"
                    self[current_option] = True # default argument w/o value
                    counter = 0        # number of values for this option
                else:                  # this entry is value for latest option
    
                    if counter < 1:
                        self[current_option] = cmd
    
        # in case, several values follow after a "-xxx" option convert dict
        # entry into list and add all elements (until the next "-") to it
                    else:
                        if counter == 1:   # there is already a value assigned
        # convert to list
                            self[current_option] = [ self[current_option] ]
        # add value to list
                        self[current_option] = self[current_option] + [cmd]
    
                    counter = counter + 1
    
        except (KeyError, UnboundLocalError), why:
            raise KeyError, "Can't resolve command line options.\n \tError:"\
                  +str(why)
    
    
    def isEmpty( self ):
        """@return: True if external commandline is empty"""
        return len( sys.argv ) == 1
        

    def parse( self, argv=sys.argv ):

        self.__parse( argv )

        if 'x' in self:
            self.parseFromFile( self['x'] )
    
        ## fill in missing default values
        self.weakupdate( self.defaults )
        
        self.__typecast()



## Test code ##

print sys.argv
s = '-i 0 0 -cut 25 -v 1.2 -o out.txt -a -b'
sys.argv += s.split()
print sys.argv

o = Commandline( cut=int, a=False, i=[1,2] )
o.parse()
