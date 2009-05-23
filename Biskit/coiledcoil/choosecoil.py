


    def parseData(self,dat = ""):
        
        try:
            lineas = open(dat,"r").readlines()
        except IOError, msg:
            raise BiskitError('cannot open score file %s.\n Error: %s' \
                              % (db, msg ) )
        
        lineas = [ l.strip() for l in lineas ]
        
        self.data = {}
        
        for l in lineas:
            
            aux = l.split()
            if aux!=[]:
                self.data[aux[0]]= {}
                
                
                for w in aux[1:]:
                    aux2 = w.split(":")
                    if aux2[0]!="Paper":
                        self.data[aux[0]][aux2[0]]= aux2[1]
                    else:
                        self.data[aux[0]][aux2[0]]= (aux2[1],aux2[2])
            
    def sameRegister(self,hept1 = "",hept2 ="",chain=""):
        """
        This function returns True when two heptads (hept1,hept2) have
        the same registers in chain (so their registers have to be 'abcdefg'
        for both.
        
        @param hept1: First heptad to compare.
        @type hept1: string
        @param hept2: Second heptad to compare.
        @type hept2: string
        @param chain: String to which this heptads belong.
        @type chain: string
        
        @return: True if both share registers.
        @rtype: boolean
        """
        a = chain.find(hept1)
        b = chain.find(hept2)
        
        return (a-b)%7 == 0 
    
    