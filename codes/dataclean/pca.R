library(data.table)
library(irlba)

args = commandArgs(trailingOnly=TRUE)
print(args)
if(length(args)!=6){
stop("Provide tissue name, AA counts, EA counts, CEU counts, CHB counts, YRI counts")
}
tissue.name = args[1]
n.aa = as.numeric(args[2])
n.ea = as.numeric(args[3])
n.ceu = as.numeric(args[4])
n.chb = as.numeric(args[5])
n.yri = as.numeric(args[6])
print(n.aa); print(n.ea)
nn=n.aa+n.ea+n.ceu/2+n.chb/2+n.yri/2
print(nn)

onelargefile = matrix(0, nrow=nn, ncol=1)

for (i in 1:22){
	print(paste(i, 'start reading'))
	#before pruning
	setwd(paste0('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/',tissue.name,'/runlamp/bychr/chr',i))
	#after pruning
	aa = setDF(fread('AAgeno.shortandfat.final'))
	ea = setDF(fread('EAgeno.shortandfat.final'))
	yri = setDF(fread('YRIgeno.shortandfat.final'))
	ceu = setDF(fread('CEUgeno.shortandfat.final'))	
	chb = setDF(fread('CHBgeno.shortandfat.final'))
	yri = sapply(yri, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
	ceu = sapply(ceu, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
	chb = sapply(chb, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
	dim(yri); dim(ceu); dim(chb)
	ind.yri = matrix(seq(1, (n.yri-1), by=2), ncol = 1)
	ind.ceu = matrix(seq(1, (n.ceu-1), by=2), ncol = 1)
	ind.chb = matrix(seq(1, (n.chb-1), by=2), ncol = 1)
	yri = matrix(t(apply(ind.yri, 2, function(t) yri[t, ] + yri[t+1, ])), nrow=(n.yri/2))
	ceu = matrix(t(apply(ind.ceu, 2, function(t) ceu[t, ] + ceu[t+1, ])), nrow=(n.ceu/2))
	chb = matrix(t(apply(ind.chb, 2, function(t) chb[t, ] + chb[t+1, ])), nrow = (n.chb/2))
	print(dim(yri))
	print(dim(ceu))
	print(dim(chb))
	x = rbind(aa, ea, yri, ceu, chb)	
	onelargefile = cbind(onelargefile, x)
}
C1 = onelargefile[1:(nn-n.chb), ]
C2 = onelargefile[1:nn, ]
print(dim(C1))
print(dim(C2))
print(dim(onelargefile))

print('create index')
index1 = c(rep('AA',n.aa), rep('EA',n.ea),rep('YRI',n.yri/2),rep('CEU',n.ceu/2))
index2 = c(rep('AA', n.aa), rep('EA', n.ea), rep('YRI', n.yri/2), rep('CEU', n.ceu/2), rep('CHB', n.chb/2))
index1 = as.factor(index1)
index2 = as.factor(index2)

setwd(paste0('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/',tissue.name,'/runlamp/pca/'))

print(dim(C1))
print(dim(C2))
print('centering..')
C1 = scale(C1)
C2 = scale(C2)

print('remove NaNs..')
bad1 = apply(C1, 2, function(x) all(is.nan(x)))
bad2 = apply(C2, 2, function(x) all(is.nan(x)))

print('bad columns with zero variance found')
C1 = C1[, !bad1]
C2 = C2[, !bad2]

print('the new dimensions removing NaNs')
print(dim(C1))
print(dim(C2))

#svd
print('svd-ing..')
s1 = irlba(C1, nu = 2, nv = 2)
s2 = irlba(C2, nu = 2, nv = 2)

pcvectors1 = cbind(C1%*%s1$v[,1], C1%*%s1$v[,2], index1[!bad1])
write.table(pcvectors1, 'pc_vectors1_beforeprune', row.names = FALSE, col.names = FALSE, quote = FALSE)
pcvectors2 = cbind(C2%*%s2$v[,1], C2%*%s2$v[,2], index2[!bad2])
write.table(pcvectors2, 'pc_vectors2_beforeprune', row.names = FALSE, col.names = FALSE, quote = FALSE)




png('pca1.png')

plot(pcvectors1[,1], pcvectors1[,2], col=c('red','blue','orange','green')[index1])
dev.off()
png('pca2.png')
plot(pcvectors2[,1], pcvectors2[,2], col=c('red','blue','purple','orange','green')[index2])
dev.off()
