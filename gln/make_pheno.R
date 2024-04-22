library(data.table)

pheno_type = 'LSimp'
path = '/home/panwei/he000176/jinwen/data/pheno/'

if(pheno_type=='LSimp'){
    pheno = fread(paste0(path,'LSimp_base_25000.csv'))
}
if(pheno_type=='prs'){
    pheno = fread(paste0(path,'prs_pheno_on_analysis_gln.txt'))
}

id_p = c('0.1','0.2','0.3','0.4','0.5','0.5B','0.6','0.7','0.8','0.9')
for(i in 1:10){
    cat('proportion:',id_p[i],'\n')
    sample_ids = fread(paste0(path,'training_data_for_base_models/ID_',id_p[i],'.txt'))
    sample_ids = sample_ids[,1,drop=F]
    colnames(sample_ids) = 'ID'
    tmp_pheno = merge(sample_ids, pheno, by = 'ID')
    
    ####### the code below checks for the correlation between the observed trait and the LS or PRS imputed trait ##########
    real = fread(paste0(path,'training_data_for_base_models/pheno_real_',id_p[i],'.csv'))
    ##### first check that the real pheno has the correct sample ids
    cat('real pheno has correct sample ids? ',all(real$ID==sample_ids$ID),'\n')
    tmp_check = merge(tmp_pheno, real, by='ID')
    cat('tmp_check dim:',dim(tmp_check),'\n')
    cat('correlation btw real pheno and imputed pheno:',cor(tmp_check$pheno.x,tmp_check$pheno.y),'\n')
    
    fwrite(tmp_pheno, paste0(path,'training_data_for_base_models/pheno_',pheno_type,'_',id_p[i],'.csv'))
}

