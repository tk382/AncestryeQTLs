#-----------------------preliminary----------------------------#
cmdArgs = commandArgs(trailingOnly = TRUE)
if(length(cmdArgs)!=7){stop('provide tissue name, chr, skip, step, block index, n.aa, n.ea')}
tissue = cmdArgs[1]
chr = as.numeric(cmdArgs[2])
skipnum = as.numeric(cmdArgs[3])
step = as.numeric(cmdArgs[4])
nn = as.numeric(cmdArgs[5])
n.aa = as.numeric(cmdArgs[6])
n.ea = as.numeric(cmdArgs[7])

library(data.table)
library(MASS)
library(matrixStats)
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

setDF<- function(x){
  if(!is.data.table(x)){
    stop("x must be a data.table")
  }
  setattr(x, "row.names", .set_row_names(nrow(x)))
  setattr(x, "class", "data.frame")
  setattr(x, "sorted", NULL)
  setattr(x, ".internal.sefref", NULL)
}
standardize = function(x,df=356){(x-mean(x))/(sd(x)*sqrt(df-1))}
get_p = function(ind){
  y.pi = y[ind]; ox.pi = ox[ind]; r1.pi = r1[ind]
  ssy.pi = t(t(og)*as.numeric(crossprod(og, y.pi)))
  r2.pi = scale(-ssy.pi+as.numeric(r1.pi))/sqrt(n.aa+n.ea-1)
  W.pi = LG-crossprod(t(ox.pi), crossprod(ox.pi, LG))
  ssw.pi = t(t(og)*colSums(og*W.pi))
  r3.pi = colScale(W.pi-ssw.pi)/sqrt(n.aa+n.ea-1)
  marcor.pi = colSums(r2.pi*r3.pi)*sqrt(n.aa+n.ea-1)
  return(c(max(abs(marcor.pi), na.rm = TRUE), which.max(abs(marcor.pi))))
}

setwd(paste0('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/',tissue))

sub=read.table('data/AAsubjects'); sub2=read.table('data/EAsubjects')
sub = sub$V1; sub2 = sub2$V1

tempexp = setDF(fread(paste0('data/bychr/chr',chr,'/exp'), skip = skipnum, nrow = step))
exp = tempexp[, -(1:4)]; #matrix
bigL = setDF(fread(paste0('data/bychr/chr',chr,'/L'), skip = skipnum, nrow = step)) #matrix
A = read.table('data/global', header = TRUE); sub=as.character(A$subject); A = A$g;
A = c(A, rep(0, n.ea)) #vector
phe = read.table('../../data/phenotypes_race.adjusted', header = TRUE, stringsAsFactors = FALSE)
phe = phe[c(match(sub, phe$SUBJID), match(sub2, phe$SUBJID)), ]

age = phe$AGE; gen = phe$GENDER
pos = read.table(paste0('data/bychr/chr',chr,'/pos'), header = TRUE)
cut = which(pos$pos>tempexp[1, 2]-1000000 & pos$pos < tempexp[nrow(tempexp),3]+1000000)
pos = pos[cut, ]
pos = pos$pos
if(length(cut)==0){
write.table(matrix(0, 3, 6), paste0('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/',tissue,'/interaction/chr',chr,'result',nn,'.txt'), col.names = TRUE, row.names = FALSE, quote = FALSE)
}else{
bigg = t(setDF(fread(paste0('data/bychr/chr',chr,'/geno'), skip = cut[1], nrow = length(cut))))
bigg = bigg[-(1:2), ]
bigg2 = colScale(bigg)/sqrt(n.aa+n.ea-1)
age = phe$AGE; gen = phe$GENDER

pc = read.table(paste0('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/',tissue,'/data/exppc'), header = FALSE)
pc1 = pc[,1]; pc2=pc[,2]
indicator = c(rep(0, n.aa), rep(1, n.ea))
gene_ids = tempexp$V4
m = length(gene_ids)
resulting.p = rep(0, m)
gamma.p = rep(0, m)
snpid = rep(0, m)
maxcor = matrix(0, m, 2)

for (ge in 1:m){
  #info
  name = gene_ids[ge]
  print(name)
  info = tempexp[ge, 1:4]; start = as.numeric(info[2]); end = as.numeric(info[3])
  #expression
  y = as.numeric(exp[ge, ])
  #position
  ind = which(pos>(start-1000000) & pos<(end+1000000))
  if(length(ind)>0){
    n = length(ind)
    newpos = pos[ind]
    #genotypes
    g = bigg2[, ind]
    #local
    L = as.numeric(bigL[ge, -(1:5)]); L = c(L, rep(0, n.ea))
    
    #get test statistic
    X = colScale(cbind(age, gen, A, L, indicator, pc1, pc2))
    ox = qr.Q(qr(X))

    tempq = g-crossprod(t(ox), crossprod(ox, g))
    r55 = apply(tempq, 2, function(x) sqrt(sum(x^2)))
    og = t(t(tempq)/r55)
    
    r1 = y - crossprod(t(ox), crossprod(ox, y))
    sy = t(og)%*%y; ssy = t(t(og)*as.numeric(sy))
    r2 = scale(-sweep(ssy,1,r1))/sqrt(n.aa+n.ea-1)
    
    LG=colScale(L*g)/sqrt(n.aa+n.ea-1)
    W = LG-ox %*% (t(ox)%*%LG)
    ssw = sweep(og, 2, colSums(og*W), FUN='*')
    r3 = colScale(W-ssw)/sqrt(n.aa+n.ea-1)
    marcor = colSums(r2*r3)*sqrt(n.aa+n.ea-1)
    maxcor[ge,] = c(max(abs(marcor), na.rm = TRUE), newpos[which.max(abs(marcor))])
    
    #permutation test
    B = 1000
    ref.cor = rep(-100, B)
    for (perm in 1:B){
      ind = sample(1:(n.aa+n.ea), replace = FALSE)
      ref.cor[perm] = get_p(ind)[1]
    }
    resulting.p[ge] = sum(ref.cor>maxcor[ge, 1])/B
    mu = mean(ref.cor); sigma = var(ref.cor); beta = mu/sigma; alpha = mu*beta;
    gamma.p[ge] = 1-pgamma(maxcor[ge, 1], shape = alpha, rate = beta)
  }else{
    maxcor[ge, ] = c(NA, NA); resulting.p[ge] = NA; gamma.p[ge] = NA
  }
}
colnames(maxcor) = c('maxcor', 'pos')
chrvec = rep(chr, length(gene_ids))
out = cbind(gene_ids, resulting.p, gamma.p, maxcor, chrvec)
colnames(out) = c('gene_ids', 'resulting.p', 'gamma.p', 'maxcor', 'pos', 'chr')
write.table(out, paste0('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/',tissue,'/interaction/chr',chr,'result',nn,'.txt'), col.names = TRUE, row.names = FALSE, quote = FALSE)


}
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

