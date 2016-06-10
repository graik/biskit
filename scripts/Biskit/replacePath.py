#!/usr/bin/env python
import sys
from Biskit import LocalPath, PDBModel, Trajectory
from Biskit.Dock import Complex, ComplexList
from Biskit.tools import load, dump, absfile

class InvalidTargetPath( Exception ):
    pass

def replacePath( path, old, new ):
    """
    analyze localpath and replace occurence of old to new.
    @param path: LocalPath, 
    @param old: str, path component to be replaced
    @param new: str, replacement
    @return: LocalPath, new LocalPath
    """
    if not path:
        return path
    
    assert isinstance( path, LocalPath )

    r = LocalPath( str(path).replace( old, new ) )

    print '\t%r -> %r' % (path, r)

    return r



def replace_in_PDBModel( o, old, new ):
    assert isinstance( o, PDBModel ), 'invalid object: %r' % o
    o.source = replacePath( o.source, old, new )
    if o.source:
        if not o.source.exists():
            raise InvalidTargetPath( 'target %r does not exist.' % o.source )

def replace_in_Complex( o, old, new ):
    assert isinstance( o, Complex ), 'invalid object: %r' % o
    replace_in_PDBModel( o.rec_model, old, new )
    replace_in_PDBModel( o.lig_model, old, new )

def replace_in_Traj( o, old, new ):
    assert isinstance( o, Trajectory ), 'invalid object: %r' % o
    replace_in_PDBModel( o.ref, old, new )

def replace_in_dict( o, old, new ):
    assert isinstance( o, dict ), 'invalid object: %r' % o
    for v in o.values():
        if isinstance( v, PDBModel ):
            replace_in_PDBModel( v, old, new )

def replace_in_list( o, old, new ):
    assert isinstance( o, list ), 'invalid object: %r' % o
    for v in o:
        if isinstance( v, PDBModel ):
            replace_in_PDBModel( v, old, new )

def replace_in_ComplexList( o, old, new ):
    assert isinstance( o, ComplexList ), 'invalid object: %r' % o
    for com in o:
        replace_in_Complex( com, old, new )

    ## The internal model cache (ComplexModelRegistry) needs to be re-created:
    o = o.__class__( list( o ) )
    return o

def repair( f, old, new ):
    """
    Load and repickle with corrected path
    """
    try:
        o = load( f )
        print 'working on %s...' % f
        if isinstance( o, PDBModel)  : replace_in_PDBModel( o, old, new )
        if isinstance( o, Complex )  : replace_in_Complex(  o, old, new )
        if isinstance( o, Trajectory): replace_in_Traj( o, old, new )
        if isinstance( o, dict )     : replace_in_dict( o, old, new )
        if isinstance( o, list )     : replace_in_list( o, old, new )
        if isinstance( o, ComplexList): o = replace_in_ComplexList(o,old,new)

        dump( o, f )
        print
        
    except InvalidTargetPath, e:
        print e
        
##########
## MAIN
##########

print """
Replace source location within pickled PDBModels, Complex, Trajectory
instances and lists or dictionaries of PDBModels.

Use::
   zsh
   cd biskit/Biskit/testdata
   replacePath.py **/*list **/*model **/*complex **/*dat **/*dic
"""

old = '/test/'
new = '/Biskit/testdata/'

pickles = sys.argv[1:]

for f in pickles:
    repair( f, old, new )



