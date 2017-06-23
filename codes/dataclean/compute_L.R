#compute local ancestry for each gene

#--------load packages--------#
library(data.table)

#-----------inputs----------#
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=3){stop('provide tissue name, AA sample size, EA sample size')}
tissue = args[1]
naa = as.numeric(args[2])
nea = as.numeric(args[3])
setwd(paste0('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/',tissue,'/data/'))

#------------read data-----------#

length_per_chr = c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468)


global = rep(0,naa)
for (ch in 1:22){
	setwd(paste0('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/',tissue,'/runlamp/bychr/chr',ch))
	averageanc = read.table('averageanc.txt')
	global = global + averageanc$V2 * length_per_chr[ch] / sum(length_per_chr)
	f = setDF(fread('finalAncestry.txt'))
	f = f[, -c(1, ncol(f))]
	pos = setDF(fread('pos.final')); pos = pos$V1;
	exp = setDF(fread(paste0('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/',tissue,'/data/bychr/chr',ch,'/exp')))
	
	exp = exp[, 1:(naa+4)]
	local = matrix(0, nrow(exp), (naa+5))
	local = as.data.frame(local); 
	colnames(local) = c('chr','start','end','gene_id','pos',colnames(exp)[5:(naa+4)])
	local[,1:3] = exp[,1:3]
	local[,4] = as.character(local[,4]); local[,4] = exp[,4]
	for (j in 1:nrow(exp)){
		where = which.min(abs(pos-(exp[j,2]+exp[j,3])/2))
		local[j,5] = pos[where]
		local[j,6:(naa+5)] = f[seq(2,(naa*2),by=2),where]*2
	}
	write.table(local, paste0('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/',tissue,'/data/bychr/chr',ch,'/L'), col.names = TRUE, row.names = FALSE, quote = FALSE)
	
}
global = data.frame(subject = colnames(exp)[5:(naa+4)], g = global)
write.table(global, paste0('/group/im-lab/nas40t2/tae/differentialeQTLs/bytissues/',tissue,'/data/global'), col.names = TRUE, row.names = FALSE, quote = FALSE)
