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
## $Revision$
## last $Date$
## last $Author$

"""
Path handling.
"""

import Biskit.tools as t
from Biskit import EHandler
from Biskit.Errors import BiskitError

import os.path
import string
import re

class LocalPathError( BiskitError ):
    pass

class LocalPath( object ):
    """
    Encapsulate a file name that might look differently in different
    environments depending on environment settings.
    The file name has an absolute (original) variant but parts of it
    can be looked up from environment variables if the original path
    doesn't exist in the current environment.

    Environment variable values must be at least 3 char long and contain
    at least one '/'. That's perhaps too restrictive? Just want to avoid
    substitution of single letters or '/a' etc.

    LocalPath tries behaving like the simple absolute path string when it
    comes to slicing etc::
      l = LocalPath( '{/home/raik|$USER}/data/x.txt' )
      l[:] == '/home/raik/data/x.txt' == l.local()
      l[-4:] == '.txt'
      str( l ) == l.local()

      l.formatted() == '{/home/raik|$USER}/data/x.txt'

    @todo: - simple-minded implementation, could be made more intelligent
    @todo: - creation from formatted string not yet implemented, e.g.::
             input:  allow multiple or overlapping substitutions
             input:  restrict allowed environment variables
             output: only substitute what is necessary until path exists
    """

    ex_fragment = re.compile('\{([a-zA-Z0-9_/ ]+)\|\$([a-zA-Z0-9_/ ]+)\}')


    def __init__( self, path=None, checkEnv=1, minLen=3, **vars ):
        """
        Create a new environment-dependent path from either a list of
        fragments and their substitution variable names or from a path or
        from a formatted string (not implemented).
        A path will be analyzed to substitute as big chunks as possible
        by environment variables. ~ and ../../ will be expanded both in
        the given path and in the environment variables.
        
        @param path: path(s)
        @type  path: [ (str, str) ] OR str
        @param checkEnv: look for substitution values among environment
                         variables (default 1)
        @type  checkEnv: 1|0
        @param vars: alternative envVar=value pairs, to be used instead of
                      environment variables
        @type  vars: envVar=value
        """

        self.fragments = [] ## list of tuples (absolut,variable_name)
        if path:
            self.set( path, checkEnv=checkEnv, minLen=minLen, **vars )

        self.__hash = None


    def __setstate__(self, state ):
        """
        called for unpickling the object.
        """
        self.__dict__ = state
        ## backwards compability
        self.__hash = getattr( self, '_LocalPath__hash', None )


    def local( self, existing=0 ):
        """
        Return a valid, absolute path. Either the existing original or with
        all substitutions for which environment variables exist.
        This function is time consuming (absfile - os.realpath is the culprit).
                         
        @param existing: don't return a non-existing path                
        @type  existing: 0|1
        
        @return: valid absolute path in current environment
        @rtype: str
        
        @raise LocalPathError: if existing==1 and no existing path can be
                               constructed via environment variables
        """
        result = string.join( [ f[0] for f in self.fragments ], '' )
        result = t.absfile( result )

        if os.path.exists( result ):
            return result

        result = ''
        for abs, env  in self.fragments:
            if env:
                result += os.environ.get( env, abs )
            else:
                result += abs

        result = t.absfile( result )

        if existing and not os.path.exists( result ):
            raise LocalPathError, "Can't construct existing path from %s."%\
                  self.formatted()

        return result


    def formatted( self ):
        """
        Get a string representation that describes the original path and all
        possible substitutions by environment variables.
        
        @return: formated path e.g. C{ '{/home/raik|$USER}/data/x.txt' }
        @rtype: str
        """
        r = ""
        for absolut,var in self.fragments:
            if var:
                r += "{%s|$%s}" % (absolut, var)
            else:
                r += absolut
        return r


    def original( self ):
        """
        Get the original path (also non-absolute) that is used if there are
        no environment variables for any of the substitutions.
        
        @return: original path
        @rtype: str
        """
        result = [ f[0] for f in self.fragments ]
        return string.join( result, '' )


    def set( self, v, checkEnv=1, minLen=3, **vars ):
        """
        Assign a new file name
        
        @param v: fragment tuples or path
        @type  v: [ (str,str) ] OR str
        @param checkEnv: look for possible substitutions in environment
        @type  checkEnv: 0|1
        @param minLen: mininal length of environment variables to consider
        @type  minLen: int
        """
        if type( v ) == list:
            self.set_fragments( v )
        else:

            if type( v ) == str and (checkEnv or vars):
                self.set_path( v, minLen=minLen, **vars )
            else:

                if type( v ) == str:
                    self.set_string( v )
                else:
                    raise PathError, 'incompatible value for PathLink' + str(v)

        self.__hash = None


    def set_fragments( self, *fragments ):
        """
        Set new path from list of path fragments and their possible
        environment variable substitutions. Fragments that can not be
        substituted are given as (str, None).
        
        @param fragments: list of fragment tuples
        @type  fragments: [ (str, str), (str, None), .. ]
        """
        self.fragments = fragments
        self.__hash = None


    def set_string( self, s ):
        """
        Set a new path and its substitutable parts from a formatted string.
        
        @param s: formatted like {/home/raik|$USER}/data/test.txt
        @type  s: str

        @raise PathError: formatted input is not yet implemented
        """
        fragments = []
        pos = 0

        for m in self.ex_fragment.finditer( s ):

            if m.start() > pos:
                fragments += [ ( s[pos:m.start()], None ) ]

            fragments += [ m.groups() ]

            pos = m.end()

        if len( s ) > pos:
            fragments += [ ( s[pos:len(s)], None ) ]

        self.set_fragments( *fragments )


    def set_path( self, fname, minLen=3, **vars ):
        """
        Set a new path and try to identify environment variables that could
        substitute parts of it. If vars is given, env. variables are ignored.
        
        @param fname: relative or absolute file name
        @type  fname: str
        @param vars: alternative param=value pairs with suggested substitutors
        @type  vars: param=value
        """
        env_items = self.__paths_in_env( minLen=minLen, vars=vars )

        fragments = [ ( fname, None ) ] ## default result

        for name, value in env_items:
            fragments = self.__substitute( fragments, name, value )

        self.fragments = fragments

        self.__hash = None


    def exists( self ):
        """
        Check if path exists
        
        @return: 1 if if current path exists
        @rtype: 1|0
        """
        return os.path.exists( self.local() )


    def load( self ):
        """
        Try to unpickle an object from the currently valid path.
        
        @return: unpickled object 
        @rtype: any
        
        @raise IOError: if file can not be found
        """
        try:
            return t.Load( self.local( existing=1 ) )
        except LocalPathError, why:
            raise IOError, "Cannot find file %s (constructed from %s)" %\
                  self.local(), str( self )


    def dump( self, o ):
        """
        Try to pickle an object to the currently valid path.
        
        @return: the absolute path to which o was pickled
        @rtype: str
        """
        try:
            f = self.local()
            t.Dump( f, o )
            return f
        except:
            t.errWriteln("Couldn't dump to %s (constructed from %s)" %\
                         self.formatted(), self.local() )
            raise


    def __substitute( self, fragments, name, value ):
        """
        Look in all not yet substituted fragments for parts that can be
        substituted by value and, if successful, create a new fragment
        
        @param fragments: fragment tuples
        @type  fragments: [ (str, str) ]
        @param name: substitution variable name
        @type  name: str
        @param value: susbtitution value in current environment
        @type  value: str
        
        @return: fragment tuples
        @rtype: [ (str, str) ]
        """
        result = []

        try:
            for abs, subst in fragments:

                if not subst:   ## unsubstituted fragment

                    v = t.absfile( value )
                    a = t.absfile( abs )

                    pos = a.find( v )

                    if pos != -1:
                        end = pos + len( v )

                        f1, f2, f3 = a[0:pos], a[pos:end], a[end:]

                        if f1:
                            result += [ (f1, None) ] ## unsubstituted head
                        result += [ (f2, name) ]     ## new substitution
                        if f3:
                            result += [ (f3, None) ] ## unsubstituted tail

                    else:
                        result += [ (abs, subst) ]
                else:
                    result += [ (abs, subst ) ]
        except OSError, why:
            EHandler.fatal("Substituting path fragments: \n" +
                                 str( fragments ) + '\nname: ' + str( name ) +
                                 '\nvalue:' + str( value ) )

        return result


    def __is_path( self, o, minLen=3 ):
        """
        Check whether an object is a path string (existing or not).
        
        @param minLen: minimal length of string o to be counted as path
        @type  minLen: int
        
        @return: 1|0
        @rtype: int
        """
        r = ( type( o ) == str and o.find('/') != -1 and len(o) >= minLen )
        if r:
            try:
                s = t.absfile( o )
                return 1
            except:
                return 0
        return 0


    def __paths_in_env( self, minLen=3, vars={} ):
        """
        Get all environment variables with at least one '/' sorted by length.
        
        @param vars: alternative param=value pairs to consider
                     instead of environment
        @type  vars: param=value
        
        @return: [ (variable name, value) ] sorted by length of value
        @rtype: [ (str,str) ]
        """
        vars = vars or os.environ

        ## all environment values with at least one '/' sorted by length
        items = vars.items()

        pairs = [ (len(v[1]), v) for v in items
                  if self.__is_path(v[1]) and v[0]!= 'PWD' and v[0]!= 'OLDPWD']

        pairs.sort()
        pairs.reverse()

        return [ x[1] for x in pairs ]


    def __str__( self ):
        """
        @return: Same as local(). string representation for print and str()
        @rtype: str
        """
        return self.local()


    def __repr__( self ):
        """
        @return: formatted output (Python representation)
        @rtype: str
        """
        return "LocalPath[ %s ]" % self.formatted()


    def __len__( self ):
        """
        Time costly when repeated many times.
        
        @return: length of file name in current environment
        @rtype: int
        """
        return len( self.local() )


    def __getslice__( self, a, b ):
        return self.local()[a:b]


    def __getitem__( self, i ):
        return self.local()[i]


    def __eq__( self, other ):
        """
        supports this == other -> 0|1
        """
        if not isinstance( other, LocalPath ):
            return 0
        return self.fragments == other.fragments


    def __ne__( self, other ):
        """
        supports this != other -> 0|1
        """
        if not isinstance( other, LocalPath ):
            return 1
        return self.fragments != other.fragments


    def __hash__( self ):
        """
        if __eq__ or __cmp__ are defined hash has to be defined too, otherwise
        the objects cannot be used as keys in dictionaries (needed for Complex-
        ModelRegistry).
        
        @return: int
        @rtype: 
        """
        if self.__hash == None:
            self.__hash = self.formatted().__hash__()

        return self.__hash



#############
##  TESTING        
#############
        
class Test:
    """
    Test class
    """
    
    def run( self, local=0 ):
        """
        run function test

        @param local: transfer local variables to global and perform
                      other tasks only when run locally
        @type  local: 1|0
        
        @return: 1
        @rtype: int
        """
        os.environ['PRJ_INTERFACES'] = '~raik/data/tb/interfaces'

        path = []

        l = LocalPath()

        ## Example 1
        l.set_fragments(
            ('/home/Bis/johan/data/tb/interfaces','PRJ_INTERFACES'),
            ('/c11/com_wet/ref.com', None) )
        path += [ 'Example 1:\n %s : %s \n'%(l.formatted(), l.local()) ]

        ## Example 2
        l.set_path( '/home/Bis/raik/data/tb/interfaces/c11/com_wet/ref.com',
                    USER='/home/Bis/raik' )
        path +=  [ 'Example 2:\n %s : %s \n'%(l.formatted(), l.local()) ]

        ## Example 3
        l.set_path( '/home/Bis/raik/data/tb/interfaces/c11/com_wet/ref.com' )
        path += [ 'Example 3:\n %s : %s \n'%(l.formatted(), l.local()) ]

        if local:
            for p in path:
                print p
            globals().update( locals() )

        return 1



    def expected_result( self ):
        """
        Precalculated result to check for consistent performance.

        @return: 1
        @rtype:  int
        """
        return 1
    
        

if __name__ == '__main__':

    test = Test()

    assert test.run( local=1 ) == test.expected_result()

