## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2016 Raik Gruenberg & Johan Leckner
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
"""
Parse and organize commandline options into a dictionary.
"""

import sys
import Biskit.tools as T
import Biskit as B
import Biskit.test as BT

class CommandlineError( B.BiskitError ):
    pass

class CommandLine(dict):
    """ 
    Parse and organize commandline options into a dictionary. Command line
    arguments are recognized by a leading '-'. The following patterns are
    supported:

        * -arg value    ... a single value assigned to argument 'arg'

        * -arg v1 v2    ... a list of values assigned to 'arg'

        * -arg -nextarg ... 'arg' is either True (given) or False (not given)

    The argument '-x' has special meaning and can be used to read arguments
    from a file (one line per argument, ala 'arg \tvalue'). However,
    arguments given at the commandline will override arguments listed in the
    -x file.

    CommandLine can be initialized with default values and types for
    each argument and, after parsing, CommandLine will type-cast default values
    and given arguments into a dictionary representation.

    A short description can be assigned to each argument and will be reported
    along with the default value by CommandLine.reportArguments().

    Use setTest() to define arguments and values for an independent test
    run of the script -- the test must be self-sufficient and can only
    rely on data that are checked into biskit/test.

    Example:
    ========

       Let's assume your program has been called like this::

         ~> myprogram.py -cut 11 -sample 12 11 -a 5 -z 11.5 -force

       Inside myprogram.py the command line parser could be set up like this::

         o = CommandLine( a=1, b='in.txt', cut=int, sample=[ float ] )

       This tells CommandLine that the default values for -a and -b are 1 and
       'in.txt', respectively. It also tells CommandLine that any argument -cut
       should be converted to an integer and any argument(s) -sample should
       be converted to a list of floats. Arguments to -a will also be converted
       to int since the default value has this type.

       The current command line is then read in like this::

         o.parse()

       and will result in a dictionary like this::

         o
         {'cut':11, 'b':'in.txt', 'sample':[ 12.0, 11.0 ], 'a':5, 
          'z':'11.5', 'force':True }

       Note that 'b' is taken unchanged from the default values and the value 
       of 'z' remains a string because CommandLine doesn't know anything about
       this kind of argument.

    Special options
    ===============

    -help  ... list description of arguments
    -x     ... read arguments and values from file
    -test  ... prepare test

    """

    def __init__( self, **kw ):
        self.defaults = {}  #: default values for arguments
        self.types = {}  #: map arguments to target types
        self.docs = {}   #: description of arguments
        self.description = '' #: description of program itself
        self.test = {}   #: set of arguments and values to run a test case
        self.test_cleanup = [] #: files to remove after test
        self.test_ok = None #: func returns test success
        self.test_tags = ['script'] #: tags classifying the test case
        self.program = '' #: will receive name of the program
        self.required = [] #: list of required arguments (show help screen)

        self.__setDefaults( **kw )


    def __parseDefault( self, k, v ):
        """A default value can be a value, a type, or a list thereof.
        @param k: key
        @type k: str
        @param v: value
        @type v: any
        @return: type and default value (or None, False) extracted from k and v
        @rtype: (type, any)
        """

        tdefault = str #: default type
        if v is not None:
            tdefault = type(v)

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
            t, param = self.__parseDefault( k, v )

            self.types[k] = t

            if param is None:
                pass
            elif param == [None]:
                self.defaults[k] = []
            else:
                self.defaults[k] = param 

    def setTest( self, fsuccess=None, tags=[ BT.SCRIPT ], **testparam ):
        """Provide arguments for a unittest.
        @param fsuccess: function to test for test success (-> True)
        @type  fsuccess: func
        @param tags: tags to associate with this test case ([ SCRIPT ])
        Note: the tags are defined in Biskit.test
        @type  tags: [ int ]

        The remaining arg=value pairs are interpreted as test arguments
        """
        self.test = testparam
        self.test_ok = fsuccess
        self.test_tags = tags

    def setTestCleanup( self, cleanuplist ):
        """
        @param cleanup: list of files to remove after a test
        @type  cleanup: [ str ]
        """
        self.test_cleanup = cleanuplist


    def setDocs( self, **docs ):
        """
        Provide a short description of each command line argument.
        @param docs - key=value pairs (argument='description')
        @type docs  - { str:str }
        """
        self.docs = docs
        self.docs.update( test='perform self-test', help='display help screen',
                          x='get extra options from this file' )

    def setDescription( self, s ):
        """
        Provide the first part of the help screen, printed before the
        listing of argument descriptions.
        @param s: description
        @type s : str
        """
        self.description = s

    def setRequired( self, l ):
        """
        CommandLine will display the help screen if a required argument is
        missing.
        @param l: list of required arguments
        @type l: [ str ]
        """
        self.required = l

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
        t = self.types.get( k, None )
        if t is None:
            return True

        if type( t ) is list:
            t = t[0]
        assert type( t ) is type

        if type( v ) is list:
            r = True        
            for i in v:
                r = r and self.isvalid( k, i )
            return r

        return type( v ) is t

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

        @return: command dictionary::
                 ala {'pdb':['in1.pdb', 'in2.pdb'], 'psf':'in.psf'}
        @rtype: {<option> : <value>}
        """
        self.program = argv[0]
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
        """
        Fetch command line arguments and match them against the default values
        and types.
        @param argv: [default: sys.argv]
        @type  argv: dict
        """
        self.__parse( argv )

        if 'x' in self:
            self.parseFromFile( self['x'] )

        ## fill in missing default values
        self.weakupdate( self.defaults )

        self.__typecast()

        self.validate()

        self.postprocess()


    def postprocess( self, helpexit=True ):
        """
        Respond to special options or missing arguments:
        * test -- replace dictionary content by test setup
        * help -- display help screen
        @param helpexit: exit program after help screen display (True)
        @type helpexit: bool
        """
        if 'test' in self:
            self.makeTest()

        if 'help' in self or self.missingArguments():
            print self.report()
            if helpexit:
                sys.exit(0)


    REPORT_LAYOUT = '  -%-8s%-15s\t%s\n'

    def reportArguments( self ):
        """Pretty-print arguments, default values and documentation
        @return: formatted string ready for print
        @rtype: str
        """
        keys = list( set( self.defaults.keys() + self.docs.keys() +\
                          self.types.keys() ) )
        keys.sort()
        r = ''

        for k in keys:

            arg = k
            doc = self.docs.get( k, '' )

            v = self.defaults.get(k, None)
            if v is None:
                v = ''
            else:
                v = '['+str(v)+']'

            r += self.REPORT_LAYOUT % (arg, v, doc)

        return r

    def missingArguments( self ):
        """
        @return: required arguments that have not been specified
        @rtype: [ str ]
        """
        return [ a for a in self.required if not a in self ]

    def reportRequired( self ):
        """
        @return: pretty-formatted list of required (and missing) arguments
        @rtype: str
        """
        r = ''
        if self.required:
            r = 'Required arguments: '
            for a in self.required:
                r += ' -%s ' %  a

        if self.missingArguments():
            r += '\nMissing arguments : '
            for a in self.missingArguments():
                r += ' -%s ' % a

        return r

    def report( self ):
        """
        @return: formatted string with program name, description,
                 argument documentation, list of required arguments
        @rtype: str
        """
        return self.program +'\n' + \
               self.description + '\n' + \
               self.reportArguments() +'\n' + \
               self.reportRequired() + '\n'

    def makeTest( self ):
        """
        Put test arguments into the CommandLine record.
        """
        self.clear()
        self.weakupdate( self.defaults )
        self.__typecast()
        self.validate()
        self.update( self.test )
        self['test'] = True

    def testCleanup( self, tree=True ):
        """
        Remove files marked as 'cleanup' by setTest(). 
        @param tree: try also removing directories
        @type tree: bool
        """
        for f in self.test_cleanup:
            T.tryRemove( f, tree=tree )

    def testSuccess( self ):
        if not self.test_ok:
            return True
        return self.test_ok()

#############
##  TESTING        
#############

import Biskit.test as BT

class Test(BT.BiskitTest):
    """Test class """

    def prepare(self):
        s = '-i 0 0 -cut 25 -v 1.2 -out out.txt -a -b'
        sys.argv += s.split()

        if self.local:
            print '\n\tinput arguments: ', s

        self.o = CommandLine( cut=int, a=False, i=[1,2], l=[ str ], d=22,
                              odd=float )
        self.o.setDocs( cut='cutoff value for input',
                        a='activate a',
                        l='input list',
                        i='starting values')
        self.o.setTest( cut=10, a=True, extra='extraoption' )


    def test_parsing(self):
        """CommandLine.parse test"""
        self.o.parse()
        o = self.o

        if self.local:
            print '\n\toutput dictionary: ', self.o

        self.assertEqual( o['i'], [0,0] )
        self.assertEqual( o['cut'], 25 )
        self.assertEqual( o['v'], '1.2')
        self.assertEqual( o['out'], 'out.txt' )
        self.assertEqual( o['a'], True )
        self.assertEqual( o['b'], True )
        self.assertEqual( o['d'], 22 )
        self.assertEqual( o['l'], [] )


    def test_doc(self):
        """CommandLine.reportArguments test"""
        r = self.o.reportArguments()

        if self.local:
            print '\n\tArgument documentation: \n', r

        pattern = self.o.REPORT_LAYOUT % ('cut', '', 'cutoff value for input')
        self.assert_( pattern in r )


    def test_test(self):
        """CommandLine.makeTest test"""
        import copy

        d = copy.copy(sys.argv)
        d.append( '-test' )
        self.o.parse( d )
        o = self.o

        if self.local:
            print '\ntest scenario: ', self.o

        self.assertEqual( o['i'], [1,2] )
        self.assertEqual( o['cut'], 10 )
        self.assertEqual( o['a'], True )
        self.assertEqual( o['extra'], 'extraoption' )


if __name__ == '__main__':

    BT.localTest()
