import sys, re
import biskit.core.oldnumeric as oldnumeric
import numpy as N
import sys
import time


def replaceimport(l):
    if 'import' in l:
        return l.replace('import numpy.oldnumeric as oldN', 'import biskit.core.oldnumeric as N0')
    return l

def valid(method):
    if method in oldnumeric.__all__:
        return True
    return hasattr(oldnumeric, method)

def report_problem(l):
    methods = ex_method.findall(l)
    invalid = [ m for m in methods if not valid(m) ]
    for m in invalid:
        print('Error [%-15s]: ' % m, l, end=' ')

def processline(l):
    """:return (str, int) - processed input line, number of replacements"""
    
    l = replaceimport(l)
    
    if ex_oldN.findall(l):
        methods = ex_method.findall(l)
        if N.all( [ valid(m) for m in methods ] ):
            return ex_oldN.sub('N0.', l), len(methods)
        else:
            report_problem(l)
            return l, 0
    else:
        return l, 0


if __name__ == '__main__':
    ## OPTIONS
    fname = sys.argv[1]
    replacefile = len(sys.argv) > 2 and sys.argv[2] == '-replacefile'
    
    if ('alter_code.py' in fname) or ('altered.py' in fname):
        print('refusing to change %s' % fname)
        sys.exit(0)
    
    foutname = fname + '.altered.py'
    if replacefile:
        foutname = fname
    
    ## Example Input
    # fname = '/data/raik/py/biskit/Biskit/PDBModel.py'
    example = '''            return oldN.zeros( (0,3), oldN.Float32 )
    '''
    
    ## find all instances of oldN.<some method or field>
    ## capture <some method or field>
    ex_method = re.compile(r'oldN\.([a-zA-Z\_0-9]+)[\(\.\s\)\:\,]')
    ex_oldN = re.compile(r'oldN\.')

    count = 0
    result = []
    for i, l in enumerate(open(fname).readlines()):
        
        l2, replacements = processline(l)
        count += replacements
        result += [l2]
    
    if count > 0:
        print('%3i replacements in %-30s' % (count, fname), end=' ')
        
        cmt = ['## numpy-oldnumeric calls replaced by custom script; %s\n' % time.strftime("%d/%m/%Y")] 
        result = cmt + result
        
        open(foutname,'w').writelines(result)
        
        print('\t-> %s' % foutname)
