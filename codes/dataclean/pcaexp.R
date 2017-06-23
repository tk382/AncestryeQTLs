library(data.table)
library(irlba)

tissue = commandArgs(trailingOnly=TRUE)[1]
setwd(paste0('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/',tissue,'/data'))
e = setDF(fread('expression_all'))
ee = as.matrix(e[,-(1:4)])
s = irlba(ee,nv=2)
vv = s$v
write.table(vv, 'exppc', col.names = FALSE, row.names = FALSE, quote = FALSE)

