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
Defines method decorators. Decorators are wrapping functions
that can preceed the definition of methods (since Python 2.4, using the '@'
symbol) and enforce, for example, type-checking, or static calls.
Before each call, the decorator function replaces the original function by
some modified version. The details are described in the Python 2.4
documentation.
"""

import types
import threading

def accept( *required, **optional ):
    """
    Deprecated -- use type annotations and mypy instead.
    Decorator function that enforces type checking on arguments of a method.

    Example::
     @accept( int, float, opt1=int, opt2=PDBModel )
     def method( n, fraction, opt2=PDBModel(), **args ):
         ...

    The leading 'self' argument of class methods is automatically accepted and
    the parent class of the method should thus B{not} be given as first type.
    For keyword arguments, None is always an accepted value.
    """
    def wrapper( f ):

        def new_f( *args, **kwds ):

            req = required

            ## for methods of classes, expect 'self' at first position
            if f.__code__.co_varnames[0] == 'self':
                ## replaced types.InstanceType as of 
                ## http://bugs.python.org/issue8206
                req = ( getattr(types, 'InstanceType', object), ) + req

            ## check obligatory arguments
            for (type, given) in zip( req, args ):

                assert isinstance(given, type), \
                       "argument %r does not match %s" % (given, type)

            ## check optional arguments
            for (name, value) in kwds.items():

                assert name in optional, \
                       "argument %s is not allowed" % name

                if value is None: continue

                assert isinstance(value, optional[name]), \
                       "argument %s does not match %s" % (name, optional[name])

            return f( *args, **kwds )

        new_f.__name__ = f.__name__
        return new_f

    return wrapper


def synchronized( f ):
    """
    Decorator function ensuring that parallel threads call the wrapped
    class method one after the other. That means, it is guaranteed
    that this method of a given object is never executed in
    parallel. However, different instances of the same class are not
    blocked and can still call the routine in parallel.

    Example::
      @synchronized
      def open_log_file( self, fname ):
          ...
          
    @note: The decorator adds (if not already present) a RLock object
    'lock' and a Condition object 'lockMsg' to the object holding this
    method.  The wrapped method can hence call self.lockMsg.wait() or
    self.lockMsg.notify/notifyAll() directly without any need for
    acquiring a lock or creating a self.lockMsg.

    For the same reason synchronized can only be applied to methods of objects.
    """

    def lock_call_release( *arg, **kw ):

        parent = arg[0]     ## get parent object
        assert isinstance( parent, types.InstanceType ), \
               'missing self argument'

        if not hasattr( parent, 'lock' ):
            parent.lock = threading.RLock()
        if not hasattr( parent, 'lockMsg' ):
            parent.lockMsg = threading.Condition()

        parent.lock.acquire()

        try:
            result = f( *arg, **kw )
        finally:

            parent.lock.release()

        return result

    return lock_call_release


#############
##  TESTING        
#############
import biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test case"""

    def test_decorators( self ):
        """decorators test"""
        import time
        from biskit import PDBModel

        class A:

            def __init__( self, id=1 ):
                self.id = id

            @accept( int, model=PDBModel, x=str )
            def simple( self, i, model=None, **arg ):
                i += 1
                return i

            @synchronized
            def report( self, j, i=1, **arg ):
                print(self.id)


        a = A()

        t = time.time()

        result = 0
        for i in range( 10000 ):

            result += a.simple( 8, x='a' )

        
        if self.local:
            print('timing: ', time.time() - t)
            globals().update( locals() )
                              
        self.assertEqual( result, 90000 )
        

if __name__ == '__main__':

    BT.localTest()
    
