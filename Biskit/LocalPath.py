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


import Biskit.tools as t
from Biskit import EHandler
from Biskit.Errors import BiskitError

import os.path
import string

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
    comes to slicing etc:
    l = LocalPath( '{/home/raik|$USER}/data/x.txt' )
    l[:] == '/home/raik/data/x.txt' == l.local()
    l[-4:] == '.txt'
    str( l ) == l.local()

    l.formatted() == '{/home/raik|$USER}/data/x.txt'

    ToDo: - simple-minded implementation, could be made more intelligent
          - creation from formatted string not yet implemented
    e.g.   - input:  allow multiple or overlapping substitutions
           - input:  restrict allowed environment variables
           - output: only substitute what is necessary until path exists
    """

    def __init__( self, path=None, checkEnv=1, minLen=3, **vars ):
        """
        Create a new environment-dependent path from either a list of
        fragments and their substitution variable names or from a path or
        from a formatted string (not implemented).
        A path will be analyzed to substitute as big chunks as possible
        by environment variables. ~ and ../../ will be expanded both in
        the given path and in the environment variables.
        path      - [ (str, str) ] OR str
        checkEnv  - 1||0, look for substitution values among environment
                    variables (default 1)
        vars      - alternative envVar=value pairs, to be used instead of
                    environment variables
        """

        self.fragments = [] ## list of tuples (absolut,variable_name)
        if path:
            self.set( path, checkEnv=checkEnv, minLen=minLen, **vars )

        self.__hash = None


    def __setstate__(self, state ):
        """called for unpickling the object."""
        self.__dict__ = state
        ## backwards compability
        self.__hash = getattr( self, '_LocalPath__hash', None )


    def local( self, existing=0 ):
        """
        Return a valid, absolute path. Either the existing original or with
        all substitutions for which environment variables exist.
        existing   - 0||1, don't return a non-existing path
        This function is time consuming (absfile - os.realpath is the culprit).
        -> str, absolute path in current environment
        !! LocalPathError, if existing==1 and no existing path can be
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
        -> str, e.g. '{/home/raik|$USER}/data/x.txt'
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
        -> str
        """
        result = [ f[0] for f in self.fragments ]
        return string.join( result, '' )


    def set( self, v, checkEnv=1, minLen=3, **vars ):
        """
        Assign a new file name
        v        - [ (str,str) ] OR str, fragment tuples or path
        checkEnv - 0||1, look for possible substitutions in environment
        minLen   - int, mininal length of environment variables to consider 
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
        fragments - [ (str, str), (str, None), .. ], list of fragment tuples
        """
        self.fragments = fragments
        self.__hash = None


    def set_string( self, s ):
        """
        Set a new path and its substitutable parts from a formatted string.
        s - str, formatted like {/home/raik|$USER}/data/test.txt
        """
        raise PathError, 'formatted input is not yet implemented'


    def set_path( self, fname, minLen=3, **vars ):
        """
        Set a new path and try to identify environment variables that could
        substitute parts of it. If vars is given, env. variables are ignored.
        fname - str, relative or absolute file name
        vars  - alternative param=value pairs with suggested substitutors
        """
        env_items = self.__paths_in_env( minLen=minLen, vars=vars )

        fragments = [ ( fname, None ) ] ## default result

        for name, value in env_items:
            fragments = self.__substitute( fragments, name, value )

        self.fragments = fragments

        self.__hash = None


    def exists( self ):
        """
        -> 1||0, 1 if if current path exists
        """
        return os.path.exists( self.local() )


    def load( self ):
        """
        Try to unpickle an object from the currently valid path.
        -> anything
        !! IOError, if file can not be found
        """
        try:
            return t.Load( self.local( existing=1 ) )
        except LocalPathError, why:
            raise IOError, "Cannot find file %s (constructed from %s)" %\
                  self.local(), str( self )
           
    def dump( self, o ):
        """
        Try to pickle an object to the currently valid path.
        -> str, the absolute path to which o was pickled
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
        fragments  - [ (str, str) ]
        name       - str, substitution variable name
        value      - str, susbtitution value in current environment
        -> [ (str, str) ]
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
        minLen - int, minimal length of string o to be counted as path
        -> 1||0
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
        vars - alternative param=value pairs to consider instead of environment
        -> [ (str,str) ], [ (variable name, value) ] sorted by length of value
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
        -> str, Same as local(). string representation for print and str()
        """
        return self.local()
        
    def __repr__( self ):
        """
        -> str, formatted output (Python representation)
        """
        return "LocalPath[ %s ]" % self.formatted()
    
    def __len__( self ):
        """
        Time costly when repeated many times.
        -> int, length of file name in current environment
        """
        return len( self.local() )

    def __getslice__( self, a, b ):
        return self.local()[a:b]

    def __getitem__( self, i ):
        return self.local()[i]


    def __eq__( self, other ):
        """supports this == other -> 0|1"""
        if not isinstance( other, LocalPath ):
            return 0
        return self.fragments == other.fragments

    def __ne__( self, other ):
        """supports this != other -> 0|1"""
        if not isinstance( other, LocalPath ):
            return 1
        return self.fragments != other.fragments

    def __hash__( self ):
        """
        if __eq__ or __cmp__ are defined hash has to be defined too, otherwise
        the objects cannot be used as keys in dictionaries (needed for Complex-
        ModelRegistry).
        -> int
        """
        if self.__hash == None:
            self.__hash = self.formatted().__hash__()

        return self.__hash
        

if __name__ == '__main__':

    os.environ['PRJ_INTERFACES'] = '~raik/data/tb/interfaces'

    l = LocalPath()

    l.set_fragments(
        ('/home/Bis/johan/data/tb/interfaces','PRJ_INTERFACES'),
        ('/c11/com_wet/ref.com', None) )

    print l.formatted(), " : ", l.local() 

    l.set_path( '/home/Bis/raik/data/tb/interfaces/c11/com_wet/ref.com',
                USER='/home/Bis/raik' )

    print l.formatted(), " : ", l.local()

    l.set_path( '/home/Bis/raik/data/tb/interfaces/c11/com_wet/ref.com' )

    print l.formatted(), " : ", l.local()

