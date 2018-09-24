#!/usr/bin/env python3

import sys

def replaceimport(l):
    r = l
    if 'import' in l:
        r = l.replace('Biskit ', 'biskit ')
        r = r.replace('Biskit.', 'biskit.')
    return r


if __name__ == '__main__':
    
    print( sys.argv )
    ## OPTIONS
    fname = sys.argv[1]
    replacefile = len(sys.argv) > 2 and sys.argv[2] == '-replacefile'

    foutname = fname + '.altered.py'
    if replacefile:
        foutname = fname

    count = 0
    result = []
    for i, l in enumerate(open(fname).readlines()):
        
        l2 = replaceimport(l)
        if l2 != l:
            count += 1
        result += [l2]
    
    if count > 0:
        print('%3i replacements in %-30s' % (count, fname), end=' ')
    
    open(foutname,'w').writelines(result)

    print('\t-> %s' % foutname)
    
