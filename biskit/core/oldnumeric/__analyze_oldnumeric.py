import re

f = open('__listoldnumeric_uses.txt')

lines = f.readlines()

r = {}

ex = re.compile(r'oldN\.([a-zA-Z\_0-9]+)[\(\.\s\)\:\,].*')

for l in lines:
    method = None
    try:
        method = ex.findall(l)[0]
        
        if not method in r:
            r[method] = 0

        r[method] += 1
    except:
        print('could not match: ', l, end=' ')

print('RESULT')
print()
for method, count in list(r.items()):
    print("%15s: %3i" % (method, count))

print('found: ', sum([c for m,c in list(r.items())]), ' calls of ', len(r), ' different methods')
