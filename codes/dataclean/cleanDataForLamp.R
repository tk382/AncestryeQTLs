#Necessary outputs for each chromosome
#Genotype files for European American population from GTEx. short and fat. 305 by N
#Genotype files for African American population from GTEx. short and fat. 51 by N
#Allele frequency files for pure population from 1000GP. One long vector length N
#SNP position files for each chromosome. One long vector of length N
#recombination rate for given the SNP position

##################### NECESSARY INPUTS ############################
args = commandArgs(trailingOnly = TRUE)
if(length(args)!=1){stop("Provide tissue name")}
tissue.name = args[1]
setwd(paste0('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/',tissue.name))

##################### READ INPUT FILES #############################

library(data.table)

#gtex genotype data
system('which gzip')
system(paste0('gzip ../../data/genotypes/',tissue.name,'_Analysis.snps.txt.gz'))
orig.gtex.genotype = setDF(fread(paste0('gzip -dc ../../data/genotypes/',tissue.name,'_Analysis.snps.txt.gz')))

#gtex phenotype data
orig.gtex.phenotype = read.table('../../data/phenotypes_race.adjusted', header = TRUE)
                          
#get EA and AA sample names
EA.samp = orig.gtex.phenotype[orig.gtex.phenotype$RACE==3, 'SUBJID']
AA.samp = orig.gtex.phenotype[orig.gtex.phenotype$RACE==2, 'SUBJID']
tissue.samp = colnames(orig.gtex.genotype)
                          
EA.samp2 = intersect(EA.samp, tissue.samp)
AA.samp2 = intersect(AA.samp, tissue.samp)
                          
write.table(EA.samp2, 'runlamp/intermediateData/EAsubjects', row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(AA.samp2, 'runlamp/intermediateData/AAsubjects', row.names = FALSE, quote = FALSE, col.names = FALSE)
                          
#gtex SNP infos
SNPnames = as.character(orig.gtex.genotype[, 1]); li = strsplit(SNPnames, '_')
chr.gtex = sapply(li, function(x) x[1]); pos.gtex = sapply(li, function(x) x[2]);
ref.gtex = sapply(li, function(x) x[3]); alt.gtex = sapply(li, function(x) x[4]);
                          
#relevant genotypes
EAgeno = orig.gtex.genotype[, EA.samp2]
AAgeno = orig.gtex.genotype[, AA.samp2]
                          
#Save intermediate files
print('saving intermediate files for gtex..')
                          
#save entire AAgeno and EAgeno with >5% allele frequency
mean = rowMeans(AAgeno); mean2 = rowMeans(EAgeno)
#where = which(mean>0.05 & mean2 > 0.05)
#write.table(AAgeno[where, ], 'runlamp/intermediateData/AAgeno.maf.checked', row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(EAgeno[where, ], 'runlamp/intermediateData/EAgeno.maf.checked', row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(cbind(chr.gtex[where], pos.gtex[where]), 'runlamp/intermediateData/positions.maf.checked', row.names = FALSE, col.names = FALSE, quote = FALSE)
                          
print('gtex writing done')
proc.time()
                          
#1000GP data
#library('VariantAnnotation')
#index = read.csv('../../data/sample_info.csv')
#YRI.samp = as.character(index[index$Population=='YRI', 'Sample'])
#CEU.samp = as.character(index[index$Population=='CEU', 'Sample'])
#CHB.samp = as.character(index[index$Population=='CHB', 'Sample'])
                          
#read YRI
#proc.time()
#print('doing yri..')
#param.yri = ScanVcfParam(info = 'AC', samples = YRI.samp)
#vcf.yri = readVcf('../../data/raw.vcf.gz', 'b37', param = param.yri)
#geno.yri  = geno(vcf.yri)$GT
#geno.yri = as.data.frame(geno.yri, row.names = FALSE)
#geno.yri[geno.yri=='.'] = '?/?'
#x = setDT(geno.yri)[, unlist(lapply(.SD, tstrsplit, split='[/]', type.convert = TRUE), recursive = FALSE)]
#YRIgeno = setDF(x)
                          
#read CEU
#proc.time()
#print('doing ceu..')
#param.ceu = ScanVcfParam(info = 'AC', samples = CEU.samp)
#vcf.ceu = readVcf('../../data/raw.vcf.gz', 'b37', param = param.ceu)
#geno.ceu = geno(vcf.ceu)$GT
#geno.ceu = as.data.frame(geno.ceu, row.names = FALSE)
#geno.ceu[geno.ceu=='.'] = '?/?'
#x = setDT(geno.ceu)[, unlist(lapply(.SD, tstrsplit, split='[/]', type.convert = TRUE), recursive = FALSE)]
#CEUgeno = setDF(x)

#read CHB
#proc.time()
#print('doing chb..')
#param.chb = ScanVcfParam(info='AC', samples = CHB.samp)
#vcf.chb = readVcf('../../data/raw.vcf.gz', 'b37', param = param.chb)
#geno.chb = geno(vcf.chb)$GT
#geno.chb = as.data.frame(geno.chb, row.names = FALSE)
#geno.chb[geno.chb=='.'] = '?/?'
#x = setDT(geno.chb)[, unlist(lapply(.SD, tstrsplit, split='[/]', type.convert = TRUE), recursive = FALSE)]
#CHBgeno = setDF(x)
                          
#print(dim(YRIgeno))
#print(dim(CEUgeno))
#print(dim(CHBgeno))
                          
#vectors of SNP infos in 
#rr = rowRanges(vcf.yri)
#pos.vcf = start(ranges(rr))
#chr.vcf = as.integer(seqnames(rr))
                          
#save intermediate Data
#write.table(YRIgeno, 'runlamp/intermediateData/YRIgeno', row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(CEUgeno, 'runlamp/intermediateData/CEUgeno', row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(CHBgeno, 'runlamp/intermediateData/CHBgeno', row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(chr.vcf, 'runlamp/intermediateData/chr.vcf', row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(pos.vcf, 'runlamp/intermediateData/pos.vcf', row.names = FALSE, col.names = FALSE, quote = FALSE)
#print('vcf writing done')
#proc.time()

YRIgeno = setDF(fread('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/Lung/runlamp/intermediateData/YRIgeno'))
CEUgeno = setDF(fread('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/Lung/runlamp/intermediateData/CEUgeno'))
CHBgeno = setDF(fread('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/Lung/runlamp/intermediateData/CHBgeno'))
chr.vcf = setDF(fread('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/Lung/runlamp/intermediateData/chr.vcf')); chr.vcf = chr.vcf$V1
pos.vcf = setDF(fread('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/Lung/runlamp/intermediateData/pos.vcf')); pos.vcf = pos.vcf$V1
                          
##################### GET INTERSECTION POSITIONS ############################
################ SEPARATE GENOTYPES AND POSITIONS BY CHR ####################
vcf.df = data.frame(chr = chr.vcf, pos = pos.vcf)
gtex.df = data.frame(chr = chr.gtex, pos = pos.gtex)
vcf.length = rep(0, 22); gtex.length = rep(0, 22)
                          
#how many snps in each chromosome
for (i in 1:22){
  vcf.length[i] = sum(vcf.df$chr==i)
  gtex.length[i] = sum(gtex.df$chr==i)
}
                          
vcf.break = cumsum(vcf.length); gtex.break = cumsum(gtex.length)
vcf.break = c(0, vcf.break); gtex.break = c(0, gtex.break)
                          
EA.geno.list = list()
AA.geno.list = list()
YRI.geno.list = list()
CEU.geno.list = list()
CHB.geno.list = list()
pos.list = list()
                          
print('starting the loop..')
                          
for (i in 1:22){
  g.tempindex = (gtex.break[i]+1):(gtex.break[i+1])
  v.tempindex = (vcf.break[i]+1):(vcf.break[i+1])
  EAgeno.bychr = EAgeno[g.tempindex, ]
  AAgeno.bychr = AAgeno[g.tempindex, ]
  YRIgeno.bychr = YRIgeno[v.tempindex, ]
  CEUgeno.bychr = CEUgeno[v.tempindex, ]
  CHBgeno.bychr = CHBgeno[v.tempindex, ]
  vcfpos.bychr = pos.vcf[v.tempindex]
  gtexpos.bychr = pos.gtex[g.tempindex]	

  g = gtex.df[gtex.df$chr==i, 2]
  v = vcf.df[vcf.df$chr==i, 2]
  f = intersect(g, v)
  index.g = which(g %in% f)
  index.g = index.g[duplicated(g[index.g])==FALSE]
  index.v = which(v %in% f)
  finalpos = f;
  EAgen = EAgeno.bychr[index.g, ]; 
  AAgen = AAgeno.bychr[index.g, ]
  YRIgen = YRIgeno.bychr[index.v, ];
  CEUgen = CEUgeno.bychr[index.v, ]
  CHBgen = CHBgeno.bychr[index.v, ]
  
  setwd(paste0('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/',tissue.name,'/runlamp/bychr/chr',i))
                            
  #remove where pure population has more than half the values missing            
  y.count = apply(YRIgen, 1, function(x) sum(is.na(x)))
  yri.index = which(y.count>(ncol(YRIgen)/2))
  c.count = apply(CEUgen, 1, function(x) sum(is.na(x)))
  ceu.index = which(c.count>(ncol(CEUgen)/2))	
                            
  #transpose data                          
  EAgeno.shortandfat = t(round(EAgen)); AAgeno.shortandfat = t(round(AAgen))
  YRIgeno.shortandfat = t(YRIgen); CEUgeno.shortandfat = t(CEUgen)
  CHBgeno.shortandfat = t(CHBgen)
  

  #remove if the pure population has something other than two main alleles
                            
  newremove.ceu = as.numeric(which(apply(CEUgeno.shortandfat, 2, function(x) sum(!x %in% c(0, 1, NA))) != 0))
  newremove.yri = as.numeric(which(apply(YRIgeno.shortandfat, 2, function(x) sum(!x %in% c(0, 1, NA))) != 0))
  newremove.chb = as.numeric(which(apply(CHBgeno.shortandfat, 2, function(x) sum(!x %in% c(0, 1, NA))) != 0))
                            
  # remove if the allele frequency in AA and EA are less than 5%
                            
  means = colMeans(AAgeno.shortandfat)/2
  means2 = colMeans(EAgeno.shortandfat)/2
  where = which(means<0.05 | means2<0.05)
                            
  #remove where positions are duplicated
  dup = which(duplicated(f))
                            
                            
  remove = Reduce(union, list(newremove.ceu, newremove.yri, newremove.chb, where, ceu.index, yri.index, dup))
                            
  if(length(remove)>0){
    AAgeno.shortandfat = AAgeno.shortandfat[, -remove]
    EAgeno.shortandfat = EAgeno.shortandfat[, -remove]
    YRIgeno.shortandfat = YRIgeno.shortandfat[, -remove]
    CEUgeno.shortandfat = CEUgeno.shortandfat[, -remove]
    CHBgeno.shortandfat = CHBgeno.shortandfat[, -remove]
    f = f[-remove]
    summary.y = colMeans(YRIgeno.shortandfat, na.rm = TRUE)
    summary.c = colMeans(CEUgeno.shortandfat, na.rm = TRUE)
    ind = duplicated(f)
    newf = f[!ind]
    }

  write.table(EAgeno.shortandfat, 'EAgeno.shortandfat.final', row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(AAgeno.shortandfat, 'AAgeno.shortandfat.final', row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(YRIgeno.shortandfat, 'YRIgeno.shortandfat.final', row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(CEUgeno.shortandfat, 'CEUgeno.shortandfat.final', row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(CHBgeno.shortandfat, 'CHBgeno.shortandfat.final', row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(summary.y, 'YRIfrequency', row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(summary.c, 'CEUfrequency', row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(f, 'pos.final', row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
                          
proc.time()
                          
