#-----------------------preliminary----------------------------#
library(data.table)
library(MASS)
library(matrixStats)
#rs.ref = setDF(fread('gzip -dc ../../../../../../tempo/rs-position.txt.gz'))
colScale = function(x) {
  if(length(dim(x))<2){
	x = (x-mean(x))/sd(x)
  }else{
  cm = colMeans(x, na.rm = TRUE)
  csd = colSds(x, center = cm)
  x = t( (t(x) - cm) / csd )
  }
  return(x)
}
standardize = function(x,df=356){(x-mean(x))/(sd(x)*sqrt(df-1))}

get_p = function(ind,ox,y,r2){
	ox.pi = ox[ind, ]
	y.pi = y[ind]
	r1.pi = standardize(crossprod(t(ox.pi), crossprod(ox.pi,y.pi)))
	return(max(abs(as.numeric(t(r2) %*% r1.pi*sqrt(355)))))
}

setwd('/group/im-lab/nas40t2/tae/differentialeQTLs/muscle_skeletal/codes/interaction')
#setwd('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/muscle_skeletal/codes/interaction')

cmdArgs = commandArgs(trailingOnly = TRUE)
chr = cmdArgs[1]

tempexp = setDF(fread(paste0('../../data/bychr/chr',chr,'/exp')))
exp = tempexp[, -(1:4)]; #matrix
bigL = setDF(fread(paste0('../../data/bychr/chr',chr,'/L'))) #matrix
bigL = cbind(bigL, matrix(0, nrow(bigL), 301))
A = read.table('../../data/global', header = TRUE); A = A$g; A = c(A, rep(0,301))
phe = read.table('../../data/phenotypes_356indiv', header = TRUE, stringsAsFactors = FALSE)
age = phe$AGE; gen = phe$GENDER
bigg = setDF(fread(paste0('../../data/bychr/chr',chr,'/geno')))
bigg2 = colScale(bigg[,-(1:2)])/sqrt(355)
pos = bigg[,2]
pc = setDF(fread('../../data/pca_exp/pc_vectors'))

print('read data')
pc1 = pc[,1]; pc2 = pc[,2]
pc3 = pc[,3]; pc4 = pc[,4]
indicator = c(rep(0, 55), rep(1, 301))
gene_ids = tempexp$gene_id
m = length(gene_ids)
result = list()
print(m)
print('start loop')

#original : age/gen/A/L
get_result = function(ge){
  name = gene_ids[ge]
  info = tempexp[ge, 1:4]; start = as.numeric(info[2]); end = as.numeric(info[3])
  #expression
  y = as.numeric(exp[ge, ])
  #position
  ind = which(pos>(start-1000000) & pos<(end+1000000))
  if(length(ind)>0){
    n = length(ind)
    newpos = pos[ind]
    #genotypes
    g = colScale(t(bigg2[ind, ]))/sqrt(355)
    if(length(ind)==1){g = t(g)}
    #local
    L = as.numeric(bigL[ge, -(1:5)])
    #get test statistic
    X = colScale(cbind(age, gen, A, L, indicator, pc1, pc2))
    ox = qr.Q(qr(X))
    temp = crossprod(t(ox), crossprod(ox, y))
    temp2 = crossprod(t(ox), crossprod(ox, g))
    r1 = standardize(y-temp)
    r2 = colScale(g-temp2)/sqrt(355)
    out = t(r2)%*%r1*sqrt(355)
    maxcor = max(abs(out))
    maxpos = pos[which.max(abs(out))]
    tempout = matrix(0,n,2)
    #permutation test
    B = 1000
    ref.cor = rep(-1000, B)
    for (perm in 1:B){
      ind = sample(1:356, replace = FALSE)
      ref.cor[perm] = get_p(ind,ox,y,r2)
    }
    resulting.p = sum(ref.cor>maxcor)/B
    mu = mean(ref.cor); sigma = var(ref.cor); beta = mu/sigma; alpha = mu*beta;
    gamma.p = 1-pgamma(maxcor, shape = alpha, rate = beta)
    return(c(maxcor,maxpos,resulting.p,gamma.p))
  }else{
    return(c(NA, NA, NA, NA))
  }
}

result = matrix(0,m,4)
for (ge in 1:m){
  result[ge,] = get_result(ge)
}

colnames(result) = c('maxcor', 'maxpos', 'resulting.p', 'gamma.p')

write.table(result, paste0('result',chr), col.names = TRUE, row.names = FALSE, quote = FALSE)

#png(paste0(chr,'_',nn,'.png'))
# i = which.min(gamma.p)
# j = which(pos==as.numeric(out[i,5]))
# expp = as.numeric(exp[i, ])
# plot(as.numeric(as.numeric(exp[i,])) ~ bigg[,j])
# ind0 = L==0;
# ind1 = L==1; 
# ind2 = L==2;
# plot(expp[ind0] ~ bigg[ind0, j], xlim = c(-0.2, 2.2),  pch = '*',xlab = 'genotype', ylab = 'expression', main='SPATC1L, rs914245')
# points(expp[ind1]~bigg[ind1,j], col = 'green')
# points(expp[ind2]~bigg[ind2, j], col = 'red')
# 
# summary(lm(expp~A+L+indicator+age+gen+pc1+pc2+bigg[,j]+L:bigg[,j]))
# c = coef(summary(lm(expp~A+L+indicator+age+gen+pc1+pc2+bigg[,j]+L:bigg[,j])))
# abline(c[1,1], c[9,1])
# abline(c[1,1]+c[3,1], c[9,1]+c[10,1], col = 'green')
# abline(c[1,1]+2*c[3,1], c[9,1]+c[10,1]*2, col = 'red')
#dev.off()

