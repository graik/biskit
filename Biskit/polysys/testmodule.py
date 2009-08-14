#~ from undeflink import undeflink
from Biskit import *
from polyfret import *
#~ from mCitrineChromophore import *
from protein import *
from assembly import *





#~ link = UndefLink(100., 50.)
#~ f = Block ('Undefined_Link',link)

#~ entity1 = FRETProtein('mCitrineA','1huy.pdb')	
#~ f1 = Block('mCitrineA',entity1)

#~ entity2 = FRETProtein('mCitrineB','1huy.pdb')
#~ f2 = Block('mCitrineB',entity2)
#~ f2.addInterval(c)
#~ f2.addInterval(d1)
#~ pol = Assembly("myFretTestProbe")

#~ print pol

#~ pol.addBlock(f)
#~ pol.addBlock(f1)
#~ pol.addBlock(f2)
#~ pol.addBlock( Block('mCitrine',FRETProtein('mCitrine','1huy.pdb'))) #but it's better to have full access to each component

#~ print pol

#~ pol.remBlock('mCitrine')
#~ c1 = Constraint("linked", entity1,link)
#~ c2 = Constraint("linked", link,entity2)

#~ print pol

#You can still work with the features
#~ f2.name = 'mCitrineC'

#~ print pol

# And even change entities with only one step
#~ entity3 = FRETProtein('mCerulean','2Q57.pdb')
#~ f2.setEntity(entity3)

#~ print pol

#defining chromophores (as an example of a subblock)
#~ chromo = mCitrineChromophore()
#~ entity2.chromophore = chromo
#~ f3 = Block(chromo.name,chromo,entity3) # chromophore is part of entity 2
#~ pol.addBlock(f3)

#~ print pol

#~ pol.run()


# Decomposition  Test

#~ pol2 = Assembly("myTestDecomposition")

#~ p = Protein("testProtein",'2Q57_clean.pdb')

#~ pol2.addBlock(Block("testProtein",p))

#~ print pol2

#~ pol2.run()

# Basic engine Test

pol3 = Assembly("BasicEngineTest")

p1 = Protein("GFP1_en",'1HUY.pdb')
p2 = Protein("GFP2_en",'2AWT.pdb')
p3 = Protein("LeuZip_en",'1JUN.pdb')


#~ undef1 = UndefProtein("first","ASFEASSDS")
#~ undef2 = UndefProtein("second","FESDFSFG")
#~ undef3 = UndefProtein("third","ERERAFDSD")
#~ undef4 = UndefProtein("fourth","SFE")

#~ pol3.addConstraint(Constraint("linked",undef3,undef4))
#~ pol3.addConstraint(Constraint("linked",undef2,undef3))
#~ pol3.addConstraint(Constraint("linked",undef1,undef2))

pol3.addBlock(Block("GFP1_b",entity=p1))
pol3.addBlock(Block("GFP2_b",entity=p2))
pol3.addBlock(Block("LeuZip_b",entity=p3))

pol3.addConstraint(Constraint("linked",p1,p3))
#~ pol3.addConstraint(Constraint("linked",p3,p1))

pol3.addEngine(BasicEngine())

print pol3

#~ print pol3.blockdict['testProtein_undef_chain_start_A']



#~ pol3.addConstraint(Constraint("linked",pol3.blockdict['testProtein_undef_chain_start_A'].entity,undef1))



pol3.run()

#~ pol4 = Assembly("FRETProtein Test")

#~ f = FRETProtein('mCerulean')
#~ pol4.addBlock(Block("myFRETProtein",f))

#~ print pol4

