#infile  <- '/cluster/work/pausch/naveen/RECOMBINATION/GRM/male/genome.eigenvec'
infile  <- snakemake@input[["infile"]]
plotfile  <- snakemake@output[["plotfile"]]

inf  <- read.table (infile, head=F)
#inf$breed  <- substr (inf [,2],1,2)
#breeds  <- c('bv','fv')
#cols  <- c('violet', 'orange')
#inf$col  <- cols [1]
#inf [inf$breed == breeds [2], "col" ]  <- cols [2] 

pdf (plotfile,height=6, width=8)
par (mar = c (6,6,4,4),cex.lab=1.5, cex.axis=1.5, las=1)
xlabs  <- 'PC1'
ylabs  <- 'PC2'
plot (inf [,3:4],  xlab=xlabs, ylab=ylabs, pch=21, col='gray', bg='darkseagreen',cex=1.5)
#legend ('bottom', horiz=T,legend=breeds, col=cols, bty='n', cex=2, pch=19)
dev.off ()
