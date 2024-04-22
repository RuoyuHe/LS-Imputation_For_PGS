library(data.table)

p = as.integer(commandArgs(trailingOnly = TRUE)[1])
datatype = commandArgs(trailingOnly = TRUE)[2]

pheno_path = '/home/panwei/he000176/jinwen/data/pheno/'
snp_path = '/home/panwei/he000176/jinwen/PRS/'

######### load dataset #########
if(p==5 && datatype=='real'){
    print('0.5B')
    y = fread(paste0(pheno_path,'training_data_for_base_models/pheno_',datatype,'_0.',p,'B.csv'))
    cov = fread(paste0(pheno_path,'training_data_for_base_models/pheno_cov_0.',p,'B.csv'))
}else{
    y = fread(paste0(pheno_path,'training_data_for_base_models/pheno_',datatype,'_0.',p,'.csv'))
    cov = fread(paste0(pheno_path,'training_data_for_base_models/pheno_cov_0.',p,'.csv'))
}


base_id = fread(paste0(snp_path,'base_id.txt'))
colnames(base_id)='ID'
ensemble_id = fread(paste0(snp_path,'ensemble_id.txt'))
colnames(ensemble_id)='ID'
test_id = fread(paste0(snp_path,'test_id.txt'))
colnames(test_id)='ID'

load(paste0(snp_path,'base.RData'))
x_base = snp
rm(snp)
x_base = as.data.frame(cbind(base_id, x_base))
base_bim = fread(paste0(snp_path,'base_bim.bim'))
colnames(base_bim) = c('CHR','ID','POSG','POS','ALT','REF')

load(paste0(snp_path,'ensemble.RData'))
x_ensemble = snp
rm(snp)
x_ensemble = as.data.frame(cbind(ensemble_id, x_ensemble))
ensemble_bim = fread(paste0(snp_path,'ensemble_bim.bim'))
colnames(ensemble_bim) = c('CHR','ID','POSG','POS','ALT','REF')

load(paste0(snp_path,'test.RData'))
x_test = snp
rm(snp)
x_test = as.data.frame(cbind(test_id, x_test))
test_bim = fread(paste0(snp_path,'test_bim.bim'))
colnames(test_bim) = c('CHR','ID','POSG','POS','ALT','REF')

######### first we calculate GWAS on the base dataset #########
n_snps = ncol(x_base) - 1
fml = as.formula(paste0('pheno ~ geno + sex + age + ',paste0('pc',1:10,collapse=' + ')))
beta = c()
beta_se = c()
p_val = c()
cat('start calculating GWAS...\n')
for(i in 1:n_snps){
    cat('gwas snp',i,'\n')
    tmp_data = x_base[,c(1,i+1)]
    colnames(tmp_data)[2] = 'geno'
    ##### impute missing genotype #####
    if(sum(is.na(tmp_data$geno))>0){
        tmp_data$geno[is.na(tmp_data$geno)] = mean(tmp_data$geno,na.rm=T)
    }
    tmp_data = merge(y,tmp_data,by='ID')
    tmp_data = merge(tmp_data,cov,by='ID')
    m = lm(fml, data = tmp_data)
    beta = c(beta, m$coef[2])
    beta_se = c(beta_se, summary(m)$coef[2,2])
    p_val = c(p_val, summary(m)$coef[2,4])
}

gwas_results = data.frame(SNP = base_bim$ID,
                          A1 = base_bim$ALT,
                          A2 = base_bim$REF,
                          BETA = beta,
                          SE = beta_se,
                          stringsAsFactors=F
                         )
n_gwas = nrow(x_base)

fwrite(gwas_results, paste0(snp_path,'tmp_gwas/',datatype,'_gwas_0.',p,'.txt'),sep=' ')

######### compute PRS #########
prs_command_path = '/home/panwei/he000176/jinwen/'
prs_command = paste0('python3 /home/panwei/fu000217/dl/LipA/7_apply/prs-cs/PRScs/PRScs.py --ref_dir=/home/panwei/fu000217/dl/LipA/7_apply/prs-cs/ldblk_ukbb_eur',
                     ' --bim_prefix=',prs_command_path,'PRS/base_bim',
                     ' --sst_file=',prs_command_path,'PRS/tmp_gwas/',datatype,'_gwas_0.',p,'.txt',
                     ' --n_gwas=',n_gwas,
                     ' --out_dir=',prs_command_path,'PRS/tmp_prs/',datatype,'/0.',p,'/prs_imputed'
                    )
prs_message = system(prs_command,intern=TRUE)

######### combine PRS results and make predictions #########
all_effects = data.frame()
for(chr in 1:22){
    effect_file=read.table(paste0(prs_command_path,'PRS/tmp_prs/',datatype,'/0.',p,'/prs_imputed_pst_eff_a1_b0.5_phiauto_chr',chr,'.txt'))
    colnames(effect_file) = c('CHR','ID','POS','A1','A2','BETA')
    all_effects = rbind(all_effects,effect_file)
}

merged_effects = merge(base_bim,all_effects, by = 'ID',sort=F)
cat('bim ALT same as prs ALT? ',all(merged_effects$ALT==merged_effects$A1),'\n')

x_ensemble_imputed <- apply(x_ensemble, 2, function(col){
  ifelse(is.na(col), mean(col, na.rm = TRUE), col)
})
#x_ensemble_imputed = as.data.frame(x_ensemble_imputed)

x_test_imputed <- apply(x_test, 2, function(col) {
  ifelse(is.na(col), mean(col, na.rm = TRUE), col)
})
#x_test_imputed = as.data.frame(x_test_imputed)

ensemble_idx = match(merged_effects$ID, ensemble_bim$ID)
cat('any NAs in idx when mathing prs effect snps with ensemble snps?',sum(is.na(ensemble_idx)),'\n')
x_ensemble_imputed = x_ensemble_imputed[,c(1,ensemble_idx + 1)]
ensemble_pred = as.vector(x_ensemble_imputed[,-1]%*%all_effects$BETA)

test_idx = match(merged_effects$ID, test_bim$ID)
cat('any NAs in idx when mathing prs effect snps with test snps?',sum(is.na(test_idx)),'\n')
x_test_imputed = x_test_imputed[,c(1,test_idx + 1)]
test_pred = as.vector(x_test_imputed[,-1]%*%all_effects$BETA)

ensemble_results = data.frame(ID = as.character(x_ensemble_imputed[,1]), 
                              pheno = ensemble_pred,
                              stringsAsFactors=F
                             )
test_results = data.frame(ID = as.character(x_test_imputed[,1]), 
                              pheno = test_pred,
                              stringsAsFactors=F
                             )

fwrite(ensemble_results,paste0(snp_path,'results/ensemble/',datatype,'_prs_0.',p,'.csv'))
fwrite(test_results,paste0(snp_path,'results/test/',datatype,'_prs_0.',p,'.csv'))
