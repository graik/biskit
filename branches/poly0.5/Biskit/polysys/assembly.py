#~ from polyfret import *
#~ import ResiduePicker

import xplortools
from constraint import *

class Assembly :
    def __init__(self,name = "myAssembly"):
        self.blocks = []
        self.constraints = []
        self.name = name
        self.engines = []
        self.blockdict={}
        
    def addBlock(self, f):
        self.blocks.append(f)
        self.blockdict[f.name]= f
        self.blocks[-1].onInsertion(self)
        
    def remBlock(self,f):
        if isinstance(f,Block):
            self.blocks.remove(f)
        else:
            self.blocks.remove(Block(f))
        
    def addConstraint(self,c):
        self.constraints.append(c)
    
    def remConstraint(self,c):
        self.constraints.remove(c)
    
    def addEngine(self,e):
        self.engines.append(e)
    
    def remEngine(self,e):
        self.engines.remove(e)
        
    def addAssembly(self,a,addEngines= False):
        
        self.blocks.extend (a.blocks)
        self.constraints.extend (a.constraints)
        
        for sth in a.blockdict:
            self.blockdict[sth] = a.blockdict[sth]
        
        if (addEngines ):
            self.engines.extend (a.engines)
        
    
    def run(self):
        #Check blocks
        reblocks =  []
        for i in range(len(self.blocks)-1):
            reblocks.append(self.blocks[i])
            
            if reblocks.count(self.blocks[i])>1 :
                print "[WARNING Assembly run] Duplicated Block("+self.blocks[i].name+ ").Merging their intervals"
                pos = reblocks.index(self.blocks[i])
                for j in self.blocks[i].intervals:
                    reblocks[pos].addBlock(blocks[i].intervals[j])
                
        #Check constraint coherence
        for i in range (len(self.constraints)-1):
            for j in range(i+1,len(self.constraints)-1):
                clash = ConstraintClash(self.constraints[i].type,self.constraints[j].type)
                if( clash in c_clashes ) :
                    #Solve possible constraint clash
                    if self.constraints[i].type == "distance":
                        if  clash == ConstraintClash("distance","distance"):
                            pass
                            
                        elif self.constraints[i].type == "distance"  and clash == ConstraintClash("distance","glued"):
                            pass
                            
                        else:
                            print "[ERROR Constraint _build] Unhandled clash descriptor (" +clash.__str__() +")." 
                    
                    elif self.constraints[i].type == "orientation":
                        pass
                        
                    elif self.constraints[i].type == "glued":
                        pass
                        
                    else:
                        print "[ERROR Constraint _build] Unhandled constraint descriptor (" +self.constraints[i].type +")." 
        
        if( self.engines != []):
            for e in self.engines:
                e.run(self)
        
        
    def __str__(self):
        mystring = ""
        mystring+= "------------------------------------\n"
        mystring+= "Assembly: "+self.name+"\n"
        mystring+= "------------------------------------\n"
        
        if self.blocks != []:
            for i in self.blocks:
                mystring+= i.__str__()+"\n"
                mystring+= "Tree view:\n"+i.getStringTree()+"\n"
                mystring+="\n"
                
        mystring+= "\n"
        if self.constraints != []:
            for i in self.constraints:
                mystring+= i.__str__()+"\n"
        mystring+= "------------------------------------\n"
        
        return mystring
            


class Engine(object):
    def __init__(self,name = "An Engine"):
        self.name = name
    
    def run (self,myAssembly):
        print " "
        print "----------"
        print "Starting "+self.name 
        print "----------"

#for managing constraints
class BasicEngine( Engine):
    def __init__ (self):
        Engine.__init__(self,"Basic Engine")
        
        
    def run(self,myAssembly):
        Engine.run(self,myAssembly)
        
        # Separate constraints.
        lconstraints = []
        for c in myAssembly.constraints:
            if c.type == "linked":
                lconstraints.append(c)
            if c.type == "distance":
                dconstraints.append(c)
        
        # Creation of the linking tree.
        # Each row represents a structure.
        print "- Creating link structure... ",
        ctree = []
        for c in lconstraints:
            appended = False
            for crow in ctree:
                total = len(crow)
                for i in range(total):
                    if c.a in  crow[i] or c.b in crow[i] :
                        crow.append(c)
                        appended = True
            if not appended:
                ctree.append([c])
        
        print "OK"

        # Ensure that there is no connection between rows (perhaps impossible?).
        print "- Link coherence check 1 (inter chain) ... ",
        for i in range(len(ctree)):
            for c in ctree[i]:
                for c2 in ctree[i+1:]:
                    if c.a in c2 or c.b in c2:
                        #error!!
                        print " X"
                        return -1 #Cannot create treelike chains (inter row)						
        print "OK"
        
        # Check that there are not treelike chains in rows.
        print "- Link coherence check 2 (intra chain) ... ",
        for crow in ctree:
            for i in range(len(crow)):
                for c in crow[i+1:]:
                    if crow[i].triangularLink(c):
                        print " X"
                        return -2 #Cannot create treelike chains (intra row)
        print "OK"
    
        # Rows must be ordered in linking direction 
        print "- Ordering elements ... ",
        o_ctree = []
        for crow in ctree:
            order = [-1]*len(crow)
            iorder = [-1]*len(crow)
            for i in range(len(crow)):
                for j in range(len(crow)):
                    if crow[i].ilinked(crow[j]):
                        iorder[i] = j
                    if crow[i].linked(crow[j]):
                        order[i] = j
                
            # order them!
            #find the first (the one with iorder == -1)
            for i in range(len(iorder)):
                if iorder[i] == -1:
                    first = i
                    
            o_ctree.append([])
            o_ctree[-1].append(crow[first])
            next = order[first]
            while next != -1:
                o_ctree[-1].append(crow[next])
                next = order[next] 
        print "OK"
        
        # And now... link the entities
        print "- Linking elements ... "
        # Creating chain id list A - Z
        letters = range(65,min(90,len(o_ctree)+65))
        chains = range(len(o_ctree))
        #~ print len(o_ctree)
        for i in range(len(o_ctree)):
            #~ print i, chains, letters
            chains[i] = chr(letters[i])
        # Link
        i = 0
        for crow in o_ctree:
            print "\t - Creating chain "+ chains[i]+ "...",
            a = crow[0].a.run()
            for c in crow:
                a = xplortools.joinProteins(a,c.b.run())
            a['chain_id'] = a.lenAtoms() * [chains[i]]
            a.writePdb(chains[i]+".pdb")
            i = i+1
            print "OK"
            #~ print a.sequence()
        
        
        #~ for crow in  o_ctree:
            #~ for c in crow:
                #~ print c 
            #~ print "\n"

