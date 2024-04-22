library(data.table)
source('/home/panwei/he000176/jinwen/PRS/read_snp.R')

dataset = commandArgs(trailingOnly = TRUE)[1]

save_path = '/home/panwei/he000176/jinwen/PRS/'

data = read_snp(dataset)

snp = data$snp_matrix

fwrite(as.data.frame(as.character(data$ID)), paste0(save_path,dataset,'_id.txt'),col.names=F)
save(snp, file = paste0(save_path,dataset,'.RData'))
fwrite(data$total_bim, paste0(save_path,dataset,'_bim.bim'),sep='\t',col.names=F)


