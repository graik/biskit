library(lattice)
#source ('http://bioconductor.org/biocLite.R')
#biocLite('lattice')
#a = read.table('DADParry_scaled')
a = read.table('tscore_scores')
b = a[,1:length(a)]
rownames(b) = a[,1]
#b1 = b[order(rownames(b)),]

names = c('1FT1@0', '1A32@0', '1YFM@0', '256B@0', '1HIW@3', '1BPD@0', '1FOS@1', '1HLO@1', '1A36@0', '1ITF@1', '4KTQ@0', '2FHA@1', '1ZME@0', '9ICY@0', '1BGW@0', '1A02@0', '1NKN@1', '1AQT@0', '1BCF@20', '1NRE@0', '1J59@0', '1A5T@0', '1MTY@1', '3KIN@0', '1OCC@3', '1AYX@3', '4BDP@0', '1AYX@4', '1PBW@1', '1CEM@0', '1PRC@2', '1RYT@0', '1UNK@1', '1KLN@1', '1JSW@7', '1XWL@0', '2KTQ@0', '1CPY@0', '1BKD@1', '1TF4@2', '1D7M@0', '3KTQ@0', '1KTQ@1', '1AIJ@2', '1AB4@1', '1AJY@0', '1CNT@0', '1ECR@0', '1BK5@9', '1FUR@1', '1JUN@0', '1BPY@0', '2A93@0', '1IK9_cropped@0', '1RPO@0', '1R48_proposmo@0', '1DIP@0', '2Z5H@1', '1SER@0', '1A2L@0', '1PYI@0', '1A93@0', '1AEW@0', '1AT3@0', '3INK@1', '1T7P@0', '1ZYM@0', '1VSG@1', '1AFR@0', '1AX8@0', '1BPZ@0', '1RUO@0', '1HUW@1', '1CII@2', '1ECM@0', '1FGK@1', '1A92_hepadelta@1', '5EAU@0', '2PVI@0', '1BCC@1', '1BGC@0', '1AVM@0', '1RCB@1', '2BPG@0', '1DBH@0', '1EFR@1', '1AM9@2', '1DPS@13', '3BCT@1', '1TMZ@0', '2VSG@0', '1BR0@0', '1YSA@0', '1SRY@1', '1HF9_BovineATPase@0', '1ALU@1', '1G6N@0', '1RCI@1', '1AU1@3', '1E2A@1', '1VLT@2', '1BG7@0', '2SPC@1', '1BT4@0')
#colnames(b) = names
levelplot(t(as.matrix(b)))

#a = read.table('homodimeric_parallel')
#b = a[,2:8]
#rownames(b) = a[,1]
#b2 = b[order(rownames(b)),]
#colnames(b2)=c('a','b','c','d','e','f','g')
#x11()
#levelplot(t(as.matrix(b2)))
#x11()
#levelplot(t(as.matrix(b1-b2)))
#write.table(b1*b2,'homodimeric_parallel')


#source ('http://bioconductor.org/biocLite.R')
#biocLite('RColorBrewer')
#display.brewer.all()
#colors = brewer.pal(9,'YlOrRd')
#levelplot(t(as.matrix(b1)),col.regions=colors)
#biocLite('marray')
#library(marray)
#colorsFunc = colorRampPalette(colors)
#levelplot(t(as.matrix(b1)),col.regions=colorsFunc(100))
#levelplot(t(as.matrix(b1)))

