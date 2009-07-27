library(lattice)
#source ('http://bioconductor.org/biocLite.R')
#biocLite('lattice')
a = read.table('DADParry_scaled')
b = a[,2:8]
rownames(b) = a[,1]
b1 = b[order(rownames(b)),]
colnames(b1)=c('a','b','c','d','e','f','g')
levelplot(t(as.matrix(b1)))

a = read.table('homodimeric_parallel')
b = a[,2:8]
rownames(b) = a[,1]
b2 = b[order(rownames(b)),]
colnames(b2)=c('a','b','c','d','e','f','g')
x11()
levelplot(t(as.matrix(b2)))
x11()
levelplot(t(as.matrix(b1-b2)))
write.table(b1*b2,'homodimeric_parallel')


source ('http://bioconductor.org/biocLite.R')
biocLite('RColorBrewer')
display.brewer.all()
colors = brewer.pal(9,'YlOrRd')
levelplot(t(as.matrix(b1-b2)),col.regions=colors)
biocLite('marray')
library(marray)
colorsFunc = colorRampPalette(colors)
levelplot(t(as.matrix(b1-b2)),col.regions=colorsFunc(100))

