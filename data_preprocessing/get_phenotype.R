require(data.table)
library(optparse)
library(futile.logger)

option_list = list(
    make_option(c("--ukb_field_id"), type="character", default=NULL, 
              help="ukb field id for phenotype", metavar="character"),
    make_option(c("--ukb_assay_id"), type="character", default=NULL, 
              help="ukb field id for phenotype assay date", metavar="character"),
    make_option(c("--save_path"), type="character", default=NULL, 
              help="Path to save results", metavar="character"),
    make_option(c("--log_file"), type="character", default=NULL, 
              help="Path to the log file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$ukb_field_id) | is.null(opt$save_path)) {
    print_help(opt_parser)
    stop("UKB field id and save path must be supplied", call.=FALSE)
}

ukb_field_id = opt$ukb_field_id # field id for the phenotype
ukb_assay_id = opt$ukb_assay_id # field id for the assay date
save_path = opt$save_path
if(opt$log_file=="" | is.null(opt$log_file)){
    system('touch /path/to/tmp/log/file')
    log_file = '/path/to/tmp/log/file'
}else{
    log_file = opt$log_file
}

flog.appender(appender.file(log_file))
flog.threshold(INFO)

flog.info("Starting 0extractliper.R script")

cov_path = 'ukb49020.tab' # ukb phenotypes file
select=c('f.eid', ukb_field_id, paste0('f.22009.0.',1:10), 'f.22001.0.0') # f.22009.0.1 - f.22009.0.1 are genetic PCs, f.22001.0.0 is genetic sex
data=fread(cov_path,select=select)
colnames(data) = c('f.eid', 'pheno', paste0('pc',1:10), 'sex')
data1<-na.omit(data)

yob = c('f.eid', 'f.34.0.0') # year of birth
yob = fread(cov_path, select = yob, data.table = F)
colnames(yob) = c('f.eid', 'yob')
yob = na.omit(yob)

flog.info("ukb assay id: %s",ukb_assay_id) # ukb field id of the date of assay for the phenotype, used to calculate age
if (ukb_assay_id=="" | is.null(ukb_assay_id)){
    # if ukb_assay_id is not provided, use the age at recruitment
    age = fread(cov_path, select = c('f.eid', 'f.21022.0.0'), data.table = F)
    colnames(age) = c('f.eid', 'age')
    age = na.omit(age)
    idx = match(data1$f.eid, age$f.eid)
    age = age[idx,]
    data1$age = age$age
}else{
    # laze variable naming, lipA_assay is actually the assay date of the phenotype of your choosing
    lipA_assay = fread(cov_path, select = c('f.eid', ukb_assay_id), data.table = F) the ph
    colnames(lipA_assay) = c('f.eid', 'lipA_assay')
    lipA_assay = na.omit(lipA_assay)
    lipA_assay$lipA_assay = year(as.Date(lipA_assay$lipA_assay))
    age = merge(yob, lipA_assay, by = 'f.eid', all = F, sort = F)
    age$age = age$lipA_assay - age$yob
    idx = which(age$age <= 0)
    if(length(idx)>0){
      age = age[-idx,]
    }
    age = data.frame(f.eid = age$f.eid, age = age$age)
    idx = match(data1$f.eid, age$f.eid)
    age = age[idx,]
    data1$age = age$age
}



# keep only white, unrelated individuals
wb<-read.table("/panfs/jay/groups/20/panwei/shared/UKBiobankIndiv/WBA_plink_keep.txt")
kinship = fread(cov_path, select = c('f.eid','f.22021.0.0'), data.table=F)
colnames(kinship)[2] = 'kinship'
kinship = na.omit(kinship)
kinship = kinship[kinship$kinship==0,]
kinship = kinship[kinship$f.eid %in% wb$V1, 'f.eid',drop=F]

id<-as.numeric(wb[,1])
id1<-data1$f.eid
index<-intersect(id,id1)
id2<-match(index,id1)
data2<-data1[id2,]
print(dim(data2))
data2 = data2[data2$f.eid %in% kinship$f.eid,]
print(dim(data2))

pheno = data2[,c('f.eid','pheno')]
colnames(pheno) = c('ID','pheno')
cov = data2[,-2]
colnames(cov)[1] = 'ID'

flog.info(paste0('Number of NAs in data: ',sum(is.na(pheno$pheno))))

# saving the phenotype, covariates (top PCs, sex, age), and sample IDs
fwrite(pheno,paste0(save_path,'/all_pheno.csv'))
save(pheno,file=paste0(save_path,'/all_pheno.RData'))

fwrite(cov,paste0(save_path,'/all_cov.csv'))
save(cov,file=paste0(save_path,'/all_cov.RData'))

fid<-data2$f.eid
fid1<-data.frame(fid,fid)
write.table(fid1,file=paste0(save_path,'/all_ids.txt'),sep=" ", col.names=F,row.names=F)

flog.info('Saved data to path: %s',save_path)

############ sample split ##############
flog.info('Starting sample splitting')

require(data.table)
# load the data we just saved
load(paste0(save_path,'/all_pheno.RData'))
load(paste0(save_path,'/all_cov.RData'))

# make sure the sample IDs match
cat('pheno and cov have mathing IDs: ',all(pheno$ID==cov$ID),'\n')

# We first split the data into GWAS ans Analysis (see paper). train and test below actually refer to GWAS and Analysis data.
# The sample size of the Analysis data is set to 70k, rest are GWAS data.
set.seed(1)
ind_num=nrow(pheno)
test_num=70000
train_num=ind_num-test_num
train_index=sample(1:ind_num,train_num,replace=F)
test_index=(1:ind_num)[-train_index]
train_ids=pheno$ID[train_index]
test_ids=pheno$ID[test_index]

save(train_ids,file=paste0(save_path,'/gwas_ids.RData'))
save(test_ids,file=paste0(save_path,'/analysis_ids.RData'))
fwrite(data.frame(train_ids),paste0(save_path,'/gwas_ids.txt'),sep=" ", col.names=F,row.names=F)
fwrite(data.frame(test_ids),paste0(save_path,'/analysis_ids.txt'),sep=" ", col.names=F,row.names=F)
fwrite(data.frame(train_ids, train_ids, stringsAsFactors = FALSE), paste0(save_path,'/gwas_ids_forGWAS.txt'), sep=" ", col.names=F, row.names=F)
fwrite(data.frame(test_ids, test_ids, stringsAsFactors = FALSE), paste0(save_path,'/analysis_ids_forGWAS.txt'), sep=" ", col.names=F, row.names=F)

truepheno_train=pheno[train_index,]
truepheno_test=pheno[test_index,]
save(truepheno_train,file=paste0(save_path,'/truepheno_gwas.RData'))
save(truepheno_test,file=paste0(save_path,'/truepheno_analysis.RData'))
fwrite(data.frame(FID=truepheno_train$ID,IID=truepheno_train$ID,pheno=truepheno_train$pheno, stringsAsFactors = FALSE),
      paste0(save_path,'/truepheno_gwas.txt'), sep=' ')
fwrite(data.frame(FID=truepheno_test$ID,IID=truepheno_test$ID,pheno=truepheno_test$pheno, stringsAsFactors = FALSE),
      paste0(save_path,'/truepheno_analysis.txt'), sep=' ')

cov_train=cov[train_index,]
cov_test=cov[test_index,]
save(cov_train,file=paste0(save_path,'/cov_gwas.RData'))
save(cov_test,file=paste0(save_path,'/cov_analysis.RData'))

flog.info('Saved sample splitting results to %s', save_path)

############ split the Analysis data into training, ensemble, and test ##############
flog.info('Starting splitting analysis into base, ensemble, and test')

data_path = save_path
cov_save_path = save_path


library(data.table)

load(paste0(data_path,'truepheno_analysis.RData'))
analysis_pheno = truepheno_test
load(paste0(data_path,'cov_analysis.RData'))
analysis_cov = cov_test

set.seed(2)

n = nrow(analysis_pheno) 
n_folds = 10
indices = sample(n) # Shuffle indices

# Assign subsets based on specified sizes, 30k for training, 30k for ensemble learning, 10k for testing
training_indices = indices[1:30000]
ensemble_indices = indices[30001:60000]
test_indices = indices[60001:70000]

training_pheno = analysis_pheno[training_indices, ]
training_cov = analysis_cov[training_indices, ]
ensemble_pheno = analysis_pheno[ensemble_indices, ]
ensemble_cov = analysis_cov[ensemble_indices, ]
test_pheno = analysis_pheno[test_indices, ]
test_cov = analysis_cov[test_indices, ]

# save results
fwrite(training_pheno,paste0(cov_save_path,'/base_pheno_real.csv'))
fwrite(training_pheno,paste0(cov_save_path,'/base_pheno.csv'))
fwrite(data.frame(FID=training_pheno$ID,IID=training_pheno$ID,pheno=training_pheno$pheno, stringsAsFactors = FALSE),
      paste0(cov_save_path,'/base_pheno.txt'), sep=' ')
training_id = data.frame(`#FID`=training_pheno$ID,IID=training_pheno$ID)
colnames(training_id) = c('#FID', 'IID')
fwrite(training_id,paste0(cov_save_path,'/base_ID.txt'),sep=' ')
fwrite(training_cov,paste0(cov_save_path,'/base_cov.csv'))

fwrite(ensemble_pheno,paste0(cov_save_path,'/ensemble_pheno.csv'))
ensemble_id = data.frame(`#FID`=ensemble_pheno$ID,IID=ensemble_pheno$ID)
colnames(ensemble_id) = c('#FID', 'IID')
fwrite(ensemble_id,paste0(cov_save_path,'/ensemble_ID.txt'),sep=' ')
fwrite(ensemble_cov,paste0(cov_save_path,'/ensemble_cov.csv'))

fwrite(test_pheno,paste0(cov_save_path,'/test_pheno.csv'))
test_id = data.frame(`#FID`=test_pheno$ID,IID=test_pheno$ID)
colnames(test_id) = c('#FID', 'IID')
fwrite(test_id,paste0(cov_save_path,'/test_ID.txt'),sep=' ')
fwrite(test_cov,paste0(cov_save_path,'/test_cov.csv'))

flog.info('Saved base, ensemble, and test to %s', cov_save_path)

# Step 2: Splitting the training set into 10 folds; 1 fold is used for the observed trait; we impute the trait for the rest. See paper for details
fold_indices = sample(nrow(training_pheno)) 
folds = cut(fold_indices, breaks=n_folds, labels=FALSE) 
training_pheno_folds = split(training_pheno, folds)
training_cov_folds = split(training_cov, folds)

p = c(0.1,0.2,0.3,0.4,0.5)
for(i in 1:5){
    tmp_pheno = data.frame()
    tmp_cov = data.frame()
    for(j in 1:i){
        tmp_pheno = as.data.frame(rbind(tmp_pheno,training_pheno_folds[[j]]))
        tmp_cov = as.data.frame(rbind(tmp_cov,training_cov_folds[[j]]))
    }
    if (nrow(tmp_pheno) == 0 || nrow(tmp_cov) == 0) {
        stop(paste("Data for fold", p[i], "is empty. Check your data and fold indices."))
    }
    
    fwrite(tmp_pheno,paste0(cov_save_path,'/training_data_for_base_models/pheno_real_',p[i],'.csv'))
    tmp_id = data.frame(`#FID`=tmp_pheno$ID,IID=tmp_pheno$ID)
    colnames(tmp_id) = c('#FID', 'IID')
    fwrite(tmp_id,paste0(cov_save_path,'/training_data_for_base_models/ID_',p[i],'.txt'),sep=' ')
    fwrite(tmp_cov,paste0(cov_save_path,'/training_data_for_base_models/pheno_cov_',p[i],'.csv')) 
    
    ############ second part ###########
    tmp_pheno = data.frame()
    tmp_cov = data.frame()
    for(j in (i+1):n_folds){
        tmp_pheno = as.data.frame(rbind(tmp_pheno,training_pheno_folds[[j]]))
        tmp_cov = as.data.frame(rbind(tmp_cov,training_cov_folds[[j]]))
    }
    if (nrow(tmp_pheno) == 0 || nrow(tmp_cov) == 0) {
        stop(paste("Data for fold", 1 - p[i], "is empty. Check your data and fold indices."))
    }
    if(i!=5){
        fwrite(tmp_pheno,paste0(cov_save_path,'/training_data_for_base_models/pheno_real_',1-p[i],'.csv'))
        tmp_id = data.frame(`#FID`=tmp_pheno$ID,IID=tmp_pheno$ID)
        colnames(tmp_id) = c('#FID', 'IID')
        fwrite(tmp_id,paste0(cov_save_path,'/training_data_for_base_models/ID_',1-p[i],'.txt'),sep=' ')
        fwrite(tmp_cov,paste0(cov_save_path,'/training_data_for_base_models/pheno_cov_',1-p[i],'.csv'))
    }else{
        fwrite(tmp_pheno,paste0(cov_save_path,'/training_data_for_base_models/pheno_real_',1-p[i],'B.csv'))
        tmp_id = data.frame(`#FID`=tmp_pheno$ID,IID=tmp_pheno$ID)
        colnames(tmp_id) = c('#FID', 'IID')
        fwrite(tmp_id,paste0(cov_save_path,'/training_data_for_base_models/ID_',1-p[i],'B.txt'),sep=' ')
        fwrite(tmp_cov,paste0(cov_save_path,'/training_data_for_base_models/pheno_cov_',1-p[i],'B.csv'))
    }
    
}


flog.info(paste0('Saved base training data to ', cov_save_path,'/training_data_for_base_models/'))
flog.info("Completed 0extractliper.R script")