setwd('/Volumes/im-lab/nas40t2/tae/differentialeQTLs/bytissues/Adipose_Subcutaneous/interaction/')
f = list.files()
result = read.table(f[12], header = TRUE, stringsAsFactors = FALSE)
for (i in 13:length(f)){
  result = rbind(result, read.table(f[i], header = TRUE, stringsAsFactors = FALSE))
}

