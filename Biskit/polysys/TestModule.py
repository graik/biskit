
from UndefLink import UndefLink
from Biskit import *
from FRETProtein import *
from mCitrineChromophore import *
from Protein import *


#~ a = Interval(1,3)
#~ b = Interval(6,8)
c = Interval(4,5)
d1 = Interval(2,5)
#~ e = Interval(5,9)
#~ print c.checkOverlap(a,b), False
#~ print c.checkOverlap(b), False
#~ print c.checkOverlap(d1),True
#~ print c.checkOverlap(a,b,d1),True
#~ print e.checkOverlap(b),True
#~ print b.checkOverlap(e)	,True
#~ print e.checkOverlap(a),False



link = UndefLink(100., 50.)
f = Block ('Undefined_Link',link)

entity1 = FRETProtein('mCitrineA','1huy.pdb')	
f1 = Block('mCitrineA',entity1)

entity2 = FRETProtein('mCitrineB','1huy.pdb')
f2 = Block('mCitrineB',entity2)
f2.addInterval(c)
f2.addInterval(d1)
pol = Assembly("myFretTestProbe")

#~ print pol

pol.addBlock(f)
pol.addBlock(f1)
pol.addBlock(f2)
pol.addBlock( Block('mCitrine',FRETProtein('mCitrine','1huy.pdb'))) #but it's better to have full access to each component

#~ print pol

pol.remBlock('mCitrine')
c1 = Constraint("linked", entity1,link)
c2 = Constraint("linked", link,entity2)

#~ print pol

#You can still work with the features
f2.name = 'mCitrineC'

#~ print pol

# And even change entities with only one step
entity3 = FRETProtein('mCerulean','2Q57.pdb')
f2.setEntity(entity3)

#~ print pol

#defining chromophores (as an example of a subblock)
chromo = mCitrineChromophore()
entity2.chromophore = chromo
f3 = Block(chromo.name,chromo,entity3) # chromophore is part of entity 2
pol.addBlock(f3)

print pol

pol.run()


pol2 = Assembly("myTestDecomposition")

p = Protein("testProtein",'2Q57.pdb')

pol2.addBlock(Block("testProtein",p))

print pol2

pol2.run();


