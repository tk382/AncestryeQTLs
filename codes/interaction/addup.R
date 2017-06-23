library(data.table)
library(MASS)
library(matrixStats)

cmdArgs = commandArgs(trailingOnly=TRUE)

if (length(cmdArgs)!=3){stop('usage : tissue, chr, nn')}

tissue = cmdArgs[1]
chr = cmdArgs[2]
nn = cmdArgs[3]

ref = setDF(fread('gzip -dc /group/im-lab/nas40t2/tempo/rs-position.txt.gz'))
getname = function(i){
	ind = match(i, res$Pos)
	res = res[ind, 3]
	return(res)
}

setwd(paste0('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/',tissue,'/interaction/'))

out = read.table(paste0('chr',chr,'result',1,'.txt'), header = TRUE)
for (i in 2:nn){
	x = read.table(paste0('chr',chr,'result',i,'.txt'), header = TRUE)
	out = rbind(out, x)	
}

out$snpid = NA
res = ref[ref$Chr==chr, ]
out$snpid = getname(out$pos)

write.table(out, paste0('resultchr',chr,'.txt'), col.names = TRUE, row.names = FALSE, quote = FALSE)



