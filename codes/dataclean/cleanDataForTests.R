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


#tissue.name = 'Lung'
setwd(paste0('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/',tissue.name))

##################### READ INPUT FILES #############################
library(data.table)

#gtex genotype data
print('read genotype')
orig.gtex.genotype = setDF(fread(paste0('gzip -dc ../../data/genotypes/',tissue.name,'_Analysis.snps.txt.gz')))
print(colnames(orig.gtex.genotype))


#gtex phenotype data
orig.gtex.phenotype = read.table('../../data/phenotypes_race.adjusted', header = TRUE)

#gtex normalized expression data
orig.gtex.exp = read.table(paste0('../../data/v6p_fastQTL_FOR_QC_ONLY/',tissue.name,'_Analysis.v6p.FOR_QC_ONLY.normalized.expr.bed'), header = TRUE, comment='',stringsAsFactors=FALSE)
orig.gtex.exp = orig.gtex.exp[orig.gtex.exp[,4]!='X', ]
#remove chromosome X and fix positions
ref.loc = setDF(fread('../../data/gene.location.ref'))
mat = match(orig.gtex.exp[,4], ref.loc[,5])
newref = ref.loc[mat,1:3]
orig.gtex.exp[,1:3] = newref
colnames(orig.gtex.exp)[1:4] = c('chr','start','end','gene_id')
colnames(orig.gtex.exp)[-(1:4)] = gsub('[.]','-',colnames(orig.gtex.exp)[-(1:4)])

#####################SAVE SAMPLE NAMES##########################                    
EA.samp = orig.gtex.phenotype[orig.gtex.phenotype$RACE==3, 'SUBJID']
AA.samp = orig.gtex.phenotype[orig.gtex.phenotype$RACE==2, 'SUBJID']
tissue.samp = colnames(orig.gtex.genotype)
EA.samp2 = intersect(EA.samp, tissue.samp)
AA.samp2 = intersect(AA.samp, tissue.samp)
write.table(EA.samp2, paste0('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/',tissue.name,'/data/EAsubjects'), row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(AA.samp2, paste0('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/',tissue.name,'/data/AAsubjects'), row.names = FALSE, quote = FALSE, col.names = FALSE)


########################GET SNP INFOS########################
#gtex SNP infos
SNPnames = as.character(orig.gtex.genotype[, 1]); li = strsplit(SNPnames, '_')
chr.gtex = sapply(li, function(x) x[1]); pos.gtex = sapply(li, function(x) x[2]);
ref.gtex = sapply(li, function(x) x[3]); alt.gtex = sapply(li, function(x) x[4]);


print('I am done by here')
                          
#relevant genotypes
EAgeno = orig.gtex.genotype[, EA.samp2]
AAgeno = orig.gtex.genotype[, AA.samp2]

print('got EAgeno and AAgeno')

#relevant expression
EAexp = orig.gtex.exp[, EA.samp2]
AAexp = orig.gtex.exp[, AA.samp2]     
AAexp = cbind(orig.gtex.exp[, 1:4], AAexp)
write.table(cbind(AAexp, EAexp), paste0('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/',tissue.name,'/data/expression_all'), col.names = TRUE, row.names = FALSE, quote = FALSE)

#>5% allele frequency
mean = rowMeans(AAgeno); mean2 = rowMeans(EAgeno)
where = which(mean>0.05 & mean2 > 0.05)
pos.matrix = cbind(chr.gtex[where], pos.gtex[where])
AAgeno = AAgeno[where,]
EAgeno = EAgeno[where,]

AAgeno = cbind(pos.matrix, AAgeno)
colnames(AAgeno)[1:2] = c('chr','pos')

for (i in 1:22){
	setwd(paste0('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/',tissue.name,'/data/bychr/chr',i))	
	exp.count = which(AAexp[,1]==i)
	EAchrexp = EAexp[exp.count, ]
	AAchrexp = AAexp[exp.count, ]
	
	gen.count = which(pos.matrix[,1]==i)
	poschr = pos.matrix[gen.count, ]
	dup = which(duplicated(poschr[,2]))
	nodup.poschr = poschr[-dup, ]
	
	AAchrgeno = AAgeno[gen.count, ]	
	EAchrgeno = EAgeno[gen.count, ]
	AAchrgeno = AAchrgeno[-dup, ]
	EAchrgeno = EAchrgeno[-dup, ]

	write.table(cbind(AAchrexp, EAchrexp), 'exp', col.names = TRUE, row.names = FALSE, quote = FALSE)
	write.table(cbind(AAchrgeno, EAchrgeno), 'geno', col.names = TRUE, row.names = FALSE, quote = FALSE)
	write.table(AAchrgeno[,1:2], 'pos', col.names = TRUE, row.names = FALSE, quote = FALSE)	
}
