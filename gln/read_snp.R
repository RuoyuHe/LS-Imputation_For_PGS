library(data.table)
library(BEDMatrix)

read_snp <- function(dataset){
    path = '/home/panwei/he000176/jinwen/'
    if(dataset == 'base'){
        sample_ids = fread(paste0(path,'data/pheno/base_ID.txt'))$IID
    }
    if(dataset == 'ensemble'){
        sample_ids = fread(paste0(path,'data/pheno/ensemble_ID.txt'))$IID
    }
    if(dataset == 'test'){
        sample_ids = fread(paste0(path,'data/pheno/test_ID.txt'))$IID
    }
    
    gwas_bim = fread(paste0(path,'GWAS_code/bim_for_gln_001.bim'))
    colnames(gwas_bim) = c('CHR','ID','POSG','POS','ALT','REF')
    
    total_snp = data.frame()
    total_bim = data.frame()
    total_fam = data.frame()
    for(chr in 1:22){
        print(chr)
        bed_path = paste0(path,'data/bed/pruned/chr',chr)
        snp = BEDMatrix(bed_path)
        bim = fread(paste0(bed_path,'.bim'))
        colnames(bim) = c('CHR','ID','POSG','POS','ALT','REF')
        fam = fread(paste0(bed_path,'.fam'))
        
        # align sample id
        idx = which(as.character(fam$V1) %in% as.character(sample_ids))
        cat('any NAs in idx when aligning sample ids? ',sum(is.na(idx)),'\n')
        fam = fam[idx,]
        snp = snp[idx,]
        
        if(chr==1) total_fam = fam
        cat('Sample IDs in chr',chr,'the same as before? ',all(as.character(fam$V1)==as.character(total_fam$V1)),'\n')
        if (any(as.character(fam$V1)!=as.character(total_fam$V1))){
            stop(paste('Sample IDs in chr',chr,'NOT the same as before!!!'))
        }
        
        # keep snps
        idx = which(bim$ID %in% gwas_bim$ID)
        bim = bim[idx,]
        snp = snp[,idx]
        cat('All sig SNPs selected? ',all(bim$ID %in% gwas_bim$ID),'\n')
        
        total_bim = rbind(total_bim, bim)
        if(chr==1){
            total_snp = snp
        }else{
            total_snp = cbind(total_snp, snp)  
        }
        
    }
    cat('bim dim:',dim(total_bim),'\n')
    cat('fam dim:',dim(total_fam),'\n')
    cat('snp dim:',dim(total_snp),'\n')
    
    # align snps
    idx = match(gwas_bim$ID, total_bim$ID)
    if (any(is.na(idx))){
        stop(paste('Total bim has different SNPs than GWAS results!'))
    }
    
    total_bim = total_bim[idx,]
    total_snp = total_snp[,idx]
    cat('All SNPs aligned? ',all(total_bim$ID == gwas_bim$ID),'\n')
    cat('total_snp dim:',dim(total_snp),'\n')
    
    return(list(ID = total_fam$V1, snp_matrix = total_snp, snp_id = total_bim$ID))
}