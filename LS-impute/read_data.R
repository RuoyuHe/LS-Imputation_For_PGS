library(data.table)
library(BEDMatrix)

read_data <- function(impute_data = 'analysis',batch_size_idx,batch_idx){
    batch_sizes = c(7500,15000,25000,35000)
    batch_size = batch_sizes[batch_size_idx]
    
    path = '/home/panwei/he000176/jinwen/'
    if(impute_data == 'analysis'){
        sample_ids = fread(paste0(path,'data/pheno/analysis_ids.txt'),header=F)$V1
    }
    if(impute_data == 'base'){
        sample_ids = fread(paste0(path,'data/pheno/base_ID.txt'))$IID
    }
    if(impute_data == 'gwas'){
        sample_ids = fread(paste0(path,'data/pheno/gwas_ids.txt'),header=F)$V1
    }
    gwas_bim = fread(paste0(path,'GWAS_code/bim_for_prs_005.bim'))
    colnames(gwas_bim) = c('CHR','ID','POSG','POS','ALT','REF')
    gwas = fread(paste0(path,'results/GWAS/GWAS_data/gwas_byPLINK_005.txt'))
    cat('gwas results and gwas bim have the same snps? ',all(gwas$SNP==gwas_bim$ID),'\n')
    
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
        fam = fam[idx,]
        snp = snp[idx,]
        
        if(chr==1) total_fam = fam
        cat('Sample IDs in chr',chr,'the same as before? ',all(as.character(fam$V1)==as.character(total_fam$V1)),'\n')
        if (any(as.character(fam$V1)!=as.character(total_fam$V1))){
            stop(paste('Sample IDs in chr',chr,'NOT the same as before!!!'))
        }
        
        # keep snps
        idx = which(bim$ID %in% gwas$SNP)
        bim = bim[idx,]
        snp = snp[,idx]
        cat('All sig SNPs selected? ',all(bim$ID %in% gwas$SNP),'\n')
        
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
    idx = match(gwas$SNP, total_bim$ID)
    if (any(is.na(idx))){
        stop(paste('Total bim has different SNPs than GWAS results!'))
    }
    
    total_bim = total_bim[idx,]
    total_snp = total_snp[,idx]
    cat('All SNPs aligned? ',all(total_bim$ID == gwas$SNP),'\n')
    cat('total_snp dim:',dim(total_snp),'\n')
    
    start = (batch_idx-1)*batch_size + 1
    end = min(nrow(total_snp), batch_idx*batch_size)
    
    cat('start:',start,'end:',end,'\n')
    
    final_ids = total_fam$V1[start:end]
    snp_batch = total_snp[start:end,]
    
    cat('snp batch dim:',dim(snp_batch),'\n')
    
    return(list(ID = final_ids, snp_matrix = snp_batch, snp_id = total_bim$ID, gwas = gwas))
}