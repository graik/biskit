import re

f = open('__listoldnumeric_uses.txt')

lines = f.readlines()

r = {}

ex = re.compile('oldN\.([a-z\_]+)[\( \.\)\r].*')

for l in lines:
    method = None
    try:
        method = ex.findall(l)[0]
    except:
        print 'could not match: ', l,

    if method and not method in r:
        r[method] = 0

    if method:
        r[method] += 1

print 'RESULT'
print
for method, count in r.items():
    print "%15s: %3i" % (method, count)

print 'found: ', sum([c for m,c in r.items()])
