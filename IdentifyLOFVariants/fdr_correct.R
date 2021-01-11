args = commandArgs(T)
dfile = args[1]
out = args[2]

data <- read.csv(dfile,header=T,sep="\t")
varinf1 <- data[,1:8]
varinf2 <- data[,10:ncol(data)]
pval <- 1-data[,9]
#FDR <- p.adjust(pval, "bonferroni")
FDR <- p.adjust(pval, "BH")
Prob <- 1-pval
outinf <- cbind(varinf1, Prob, FDR, varinf2)
write.table(outinf,out,row.names=F,sep="\t",quote=F )
