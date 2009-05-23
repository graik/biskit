


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
            