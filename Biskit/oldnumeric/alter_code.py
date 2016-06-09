import sys, re
import Biskit.oldnumeric as oldnumeric
import numpy as N

fname = sys.argv[1]
# fname = '/data/raik/py/biskit/Biskit/PDBModel.py'
foutname = fname + '.altered.py'

example = '''            return oldN.zeros( (0,3), oldN.Float32 )
'''

## find all instances of oldN.<some method or field>
## capture <some method or field>
ex_method = re.compile('oldN\.([a-zA-Z\_0-9]+)[\(\.\s\)\:\,]')
ex_oldN = re.compile('oldN\.')

def valid(method):
    if method in oldnumeric.__all__:
        return True
    return hasattr(oldnumeric, method)

def report_problem(l):
    methods = ex_method.findall(l)
    invalid = [ m for m in methods if not valid(m) ]
    for m in invalid:
        print 'Error [%-15s]: ' % m, l,

def processline(l):
    """@return (str, int) - processed input line, number of replacements"""
    if ex_oldN.findall(l):
        methods = ex_method.findall(l)
        if N.all( [ valid(m) for m in methods ] ):
            return ex_oldN.sub('N0.', l), len(methods)
        else:
            report_problem(l)
            return l, 0
    else:
        return l, 0

count = 0
result = []
for i, l in enumerate(open(fname).readlines()):
    
    l2, replacements = processline(l)
    count += replacements
    result += [l2]

if count > 0:
    print '%3i replacements in %-30s' % (count, fname),
    open(foutname,'w').writelines(result)
    print '\t-> %s' % foutname
