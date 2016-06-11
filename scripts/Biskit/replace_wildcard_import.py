#!/usr/bin/env python
##
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
##

import Biskit.tools as t
import re, shutil

def syntax():
    print \
    """
    replace_wildcard_import -i |src_file| -m |module_name|
                            [-as |import_as| -e |exclude_method1 ..| ]

    example:
    replace_wildcard_import -i tools.py -m Numeric -as N -e Complex

    will replace all calls to any Numeric method (except Complex) by
    N.|method|
    """


def re_lst( module, exclude ):
    """
    generate list of regex that capture any call to any method of module
    """
    exec( 'methods = dir(' + module + ')' )

    r = [ re.compile( '[^a-zA-Z0-9_\.]('+m+'\()' ) for m in methods
          if not m in exclude ]

    return r


def replace_line( l, method_patterns, prefix ):
    """
    add prefix to all occurences of any method
    """
    replaced = 0
    for regex in method_patterns:
        matches = regex.findall( l )

        if matches:
            l = l.replace( matches[0], prefix + matches[0] )

            replaced += 1

    return l, replaced
    
def replace_import_statement( l, module, importas ):

    if importas == module:
        importas = ''
    else:
        importas = ' as ' + importas

    return l.replace( 'from '+module+' import *', 'import '+module+importas)


options = t.cmdDict( {'m':'Numeric', 'as':'N', 'e':['Complex'] } )

module = options['m']

importas = options.get('as', module )

srcfiles = t.toList( options.get('i', None) )

exclude = t.toList( options['e'] )


try:
    exec('import '+module)

    for fname in srcfiles:

        fname = t.absfile( fname )
        shutil.copy( fname, fname + '_' )

        methods = re_lst( module, exclude )

        fold = open( fname + '_' )
        fnew = open( fname, 'w' )

        i = 0
        for l in fold:

            i += 1

            l = replace_import_statement( l, module, importas )
            l, occurrences = replace_line( l, methods, importas + '.' )

            if occurrences > 0:
                t.errWriteln( '%s %5i %2i matches:\n\t%s' %
                              (t.stripFilename(fname), i, occurrences, l) )

            fnew.write( l )

        fnew.close()
        fold.close()

except:
    syntax()
