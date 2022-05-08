setwd("~/TWAS_testing/Rcode/power/")
library(plink2R)
library(dplyr)
library(tidyr)
library(data.table)
source("./allele_qc.R")
source("./stage2.R")

args = commandArgs(trailingOnly=TRUE)
chr = args[1]
print(paste("chr:",chr))

#GTEX
### Covariates
COV5TABLE = fread("~/deepRIV/UKB/data/GTEx/Whole_Blood.v8.covariates.txt",header = T)
### Gene Expression
gene_gtf = rtracklayer::import('/home/panwei/shared/GTEx-v8/gencode.v26.GRCh38.genes.gtf')
gene_gtf = as.data.frame(gene_gtf) %>% select(gene_id,gene_name) %>% distinct()
GENE_EXP = fread("~/deepRIV/UKB/data/GTEx/Whole_Blood.v8.normalized_expression.bed.gz")
colnames(GENE_EXP)[1] = 'chrs'
GENE_EXP$chrs = as.numeric(gsub("[^0-9.-]", "", GENE_EXP$chrs))
idx = match(GENE_EXP$gene_id, gene_gtf$gene_id)
GENE_EXP$gene_name = gene_gtf$gene_name[idx]
GENE_EXP = as.data.frame(na.omit(GENE_EXP))
gene_exp_chr = GENE_EXP %>% filter(chrs == chr)
whole_blood_buddy = '/home/panwei/he000176/deepRIV/UKB/data/Whole_Blood_buddy.txt'


### ukb phenotype
white_unrelated_keep = fread('~/deepRIV/UKB/data/white_unrelated_keep_ind.txt',
                             header = F)
ukb_pheno_all = fread('~/deepRIV/UKB/data/ukb_hdl.txt')
ukb_pheno_all = na.omit(ukb_pheno_all)
keep_idx = sort(na.omit(match(white_unrelated_keep$V1, ukb_pheno_all$f.eid)))
ukb_pheno_all = ukb_pheno_all[keep_idx,]
rm(keep_idx)
ukb_covariates = fread('~/deepRIV/UKB/data/covariates.csv')


real_data_result = list()
cat(nrow(gene_exp_chr),"\n")
for(gene_ind in 1:nrow(gene_exp_chr)){
  # get covariates
  cov5table = COV5TABLE
  # get gene expression
  gene_exp = gene_exp_chr[gene_ind,]
  
  ###UKB
  bed_dir = '/home/panwei/shared/UKBiobankIndiv/genetic_data/ukb_cal_chr'
  prefix = paste0(bed_dir,chr,'_v2')
  
  # start analysis ----------------------------------------------------------
  ukb_pheno = ukb_pheno_all
  gene_name = as.character(gene_exp$gene_id[1])
  start = floor(gene_exp$start[1]/1000) - 100
  end = ceiling(gene_exp$end[1]/1000) + 100
  
  plink_command = paste0("module load plink; \n","plink --bfile /home/panwei/he000176/deepRIV/UKB/data/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_MAF01",
                         " --chr ",chr," --from-kb ",
                         start," --to-kb ",end,
                         " --keep ",whole_blood_buddy,
                         " --geno 0 --maf 0.05 --hwe 0.001  ",
                         " --make-bed --out ~/TWAS_testing/data/For_Individual_Genes/",
                         gene_name,sep="")
  plink_msg=system(plink_command,intern=TRUE)
  snp = tryCatch(read_plink(paste("~/TWAS_testing/data/For_Individual_Genes/",gene_name,sep="")),
                 error=function(e){cat("ERROR :",
                                       conditionMessage(e),
                                       "\n")})
  if(is.null(snp))
  {
    FINAL_RESULT = NULL
    real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    
    next()
  }
  remove_command = paste("rm ~/TWAS_testing/data/For_Individual_Genes/",
                         gene_name,".*",sep="")
  system(remove_command)
  
  write.table(snp$bim$V2,paste0('~/TWAS_testing/data/ForUKB/',gene_name,'_rs.txt'),row.names=F,quote=F,col.names=F)
  
  grep_command = paste0("zgrep -F -m ",length(snp$bim$V2)," -f ~/TWAS_testing/data/ForUKB/",gene_name,"_rs.txt /home/panwei/he000176/deepRIV/UKB/data/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt  > ~/TWAS_testing/data/ForUKB/",gene_name,"_rs1.txt")
  system(grep_command)
  gtex_dic = read.table(paste0('~/TWAS_testing/data/ForUKB/',gene_name,'_rs1.txt'))
  write.table(gtex_dic$V7,paste0('~/TWAS_testing/data/ForUKB/',gene_name,'_rs.txt'),row.names=F,quote=F,col.names=F)
  
  # get UKB SNPs
  bed_dir = '/home/panwei/shared/UKBiobankIndiv/genetic_data/ukb_cal_chr'
  prefix = paste0(bed_dir,chr,'_v2')
  ukb_out = paste0('~/TWAS_testing/data/ForUKB/',gene_name)
  ukb_plink_command = paste0('module load plink; plink --bfile ',prefix,' --chr ',chr,
                             ' --keep /home/panwei/he000176/deepRIV/UKB/data/white_unrelated_keep_ind.txt',
                             ' --extract ','~/TWAS_testing/data/ForUKB/',gene_name,'_rs.txt',
                             ' --mind 0 --make-bed',
                             ' --out ',ukb_out)
  ukb_plink_msg = system(ukb_plink_command,intern=TRUE)
  
  ukb_snp = tryCatch(read_plink(paste0("~/TWAS_testing/data/ForUKB/",gene_name)),
                     error=function(e){cat("ERROR :",
                                           conditionMessage(e),
                                           "\n")})
  if(is.null(ukb_snp))
  {
    FINAL_RESULT = NULL
    real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    
    next()
  }
  remove_command = paste0("rm ~/TWAS_testing/data/ForUKB/",
                          gene_name,".*")
  system(remove_command)
  remove_command = paste0("rm ~/TWAS_testing/data/ForUKB/",
                          gene_name,"_rs*")
  system(remove_command)
  
  ### overlap between GTEx SNP and UKB
  #colnames(snp$bed) = snp$bim$V4
  rownames(ukb_snp$bed) = ukb_snp$fam[,2]
  
  #colnames(ukb_snp$bim)[4:6] = c("Position","A1","A2")
  gtex_bim_tmp = merge(snp$bim,gtex_dic,by.x='V2',by.y='V1',sort=F)
  snp$bim[,2] = gtex_bim_tmp$V7
  colnames(snp$bed) = snp$bim[,2]
  rownames(snp$bed) = snp$fam[,2]
  snp_bim_ukb = merge(snp$bim,ukb_snp$bim,by=c("V1","V2"))
  remove_flip = allele.qc(snp_bim_ukb$V5.x,snp_bim_ukb$V6.x,
                          snp_bim_ukb$V5.y,snp_bim_ukb$V6.y)
  flip_snp = snp_bim_ukb$V2[remove_flip$flip]
  snp_bim_ukb = snp_bim_ukb[remove_flip$keep,]
  
  ukb_snp_bed = ukb_snp$bed[,which(colnames(ukb_snp$bed) %in% snp_bim_ukb$V2),drop=FALSE] 
  ukb_snp_bed[,which(colnames(ukb_snp_bed) %in% flip_snp)] = 2 - ukb_snp_bed[,which(colnames(ukb_snp_bed) %in% flip_snp)]
  
  snp_bed_ind = NULL
  snp_bed_colname = colnames(snp$bed)
  for(i in 1:nrow(snp_bim_ukb))
  {
    snp_bed_ind = c(snp_bed_ind,which(snp_bed_colname == colnames(ukb_snp_bed)[i]))
  }
  SNP_BED = snp$bed[,snp_bed_ind,drop=FALSE]
  cov5table = cov5table[match(rownames(SNP_BED),cov5table$ID),]
  
  keep_ukb_indiv = intersect(ukb_pheno$f.eid,ukb_snp$fam[,2])
  ukb_pheno = ukb_pheno %>% filter(f.eid %in% keep_ukb_indiv)
  ukb_snp_bed = ukb_snp_bed[which(rownames(ukb_snp_bed)%in%keep_ukb_indiv),,drop=FALSE]
  ukb_pheno = ukb_pheno[match(rownames(ukb_snp_bed),ukb_pheno$f.eid),]
  
  if(nrow(ukb_pheno)< 3){
    FINAL_RESULT = NULL
    real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    next()
  }
  if(length(snp_bed_ind) < 2)
  {
    FINAL_RESULT = NULL
    real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    next()
  }
  ### prune
  cor_cutoff = 0.8
  cor_bed = abs(cor(SNP_BED))
  cor_bed = (cor_bed < cor_cutoff)^2
  diag(cor_bed) = 1
  i = 1
  while(i < nrow(cor_bed) )
  {
    ind = which(cor_bed[i,] == 1)
    cor_bed = as.matrix(cor_bed[ind,ind])
    i = i + 1
  }
  if(nrow(cor_bed) == 1)
  {
    FINAL_RESULT = NULL
    real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    next()
  }
  ind = which(is.element(colnames(SNP_BED),colnames(cor_bed)))
  SNP_BED = SNP_BED[,ind,drop=F]
  ukb_snp_bed = ukb_snp_bed[,ind]
  snp_bim_ukb = snp_bim_ukb[ind,,drop=F]
  
  
  ### regress gene exp on covariates
  X = as.matrix(gene_exp[1,-(1:4)])
  X = X[,match(rownames(SNP_BED),colnames(X)),drop=F]
  X1 = scale(as.numeric(X))
  
  lm_RegCov = lm(X1 ~.,data = data.frame(X1,cov5table[,-c(1)]))
  X1 = lm_RegCov$residuals
  X1 = scale(X1)
  
  if(ncol(SNP_BED) > 50)
  {
    ind = order(abs(cor(X1,SNP_BED)),decreasing = T)[1:50]
    ind = sort(ind)
    SNP_BED = SNP_BED[,ind]
    ukb_snp_bed = ukb_snp_bed[,ind]
    snp_bim_ukb = snp_bim_ukb[ind,]
  }
  
  SNP_BED = scale(SNP_BED)
  ukb_snp_bed = scale(ukb_snp_bed)
  
  ### regress ukbb pheno on covariates
  # cov_idx = match(ukb_pheno$f.eid, ukb_covariates$f.eid)
  # covariates = ukb_covariates[cov_idx,]
  # print(paste("ukb & covariates aligned:", all(covariates$f.eid==ukb_pheno$f.eid)))
  # impute missing values
  # for(i in 2:ncol(covariates)){
    # if(i==3){
      # idx = which(is.na(covariates[,i]))
      # mode_impute = as.integer(round(sum(covariates[,i], na.rm = T)/nrow(covariates)))
      # covariates[idx,i] = mode_impute
    # }else{
      # idx = which(is.na(covariates[,i]))
      # covariates[idx,i] = mean(covariates[,i], na.rm = T)
    # }
  # }
  # tmp_cov= lm(as.numeric(unlist(ukb_pheno[,2])) ~ as.matrix(covariates[,-1]))
  # ukb_pheno$no_cov = tmp_cov$residuals
  
  ### remove na
  lm_stage1_X1 = lm(X1 ~ SNP_BED)
  hatbetaX1 = lm_stage1_X1$coefficients[-1]
  na_ind = which(!is.na(hatbetaX1))
  check_ukb = apply(ukb_snp_bed,2,sd)
  na_ind_ukb = which(!is.na(check_ukb))
  na_ind = intersect(na_ind,na_ind_ukb)
  SNP_BED = SNP_BED[,na_ind,drop=FALSE]
  ukb_snp_bed = ukb_snp_bed[,na_ind,drop=FALSE]
  
  
  lm_stage1_X1 = 
    step(lm(X1~.,data = data.frame(X1,SNP_BED)),direction = "backward",trace=FALSE)
  AIC_stage1_X1 = AIC(lm_stage1_X1)
  BIC_stage1_X1 = BIC(lm_stage1_X1)
  lm_stage1_X1 = summary(lm_stage1_X1)
  stage1_sigma = lm_stage1_X1$sigma
  rsq_stage1_X1 = lm_stage1_X1$r.squared
  adjrsq_stage1_X1 = lm_stage1_X1$adj.r.squared
  se_betaX1 = hatbetaX1 = rep(0,ncol(SNP_BED))
  coef_X1 = lm_stage1_X1$coefficients
  name_coef_X1 = substr(rownames(coef_X1),1,30)
  name_SNP_BED = colnames(SNP_BED)
  for(beta_ind in 1:nrow(coef_X1))
  {
    ii = which(name_SNP_BED == name_coef_X1[beta_ind])
    hatbetaX1[ii] = coef_X1[beta_ind,1]
    se_betaX1[ii] = coef_X1[beta_ind,2]
  }
  
  num_snps = sum(hatbetaX1 != 0)
  
  if(sum(abs(hatbetaX1)) == 0 || num_snps < 2)
  {
    FINAL_RESULT = NULL
    real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    next()
  }
  
  if(nrow(ukb_snp_bed) < 500){
    FINAL_RESULT = NULL
    real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    next()
  }
  
  if(is.null(lm_stage1_X1$fstatistic)){
    FINAL_RESULT = NULL
    real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    next()
  }
  
  ### UKB
  ukb_pheno[,2] = scale(ukb_pheno[,2])
  # ukb_pheno[,3] = scale(ukb_pheno[,3])
  
  Xhat_ukb_total = ukb_snp_bed %*% hatbetaX1 + coef_X1[1,1]
  Xhat_ukb_total = scale(Xhat_ukb_total)
  y_ukb_total = as.numeric(unlist(ukb_pheno[,2]))
  # y_ukb_total = as.numeric(unlist(ukb_pheno[,3]))
  stage2_test = lm(y_ukb_total~Xhat_ukb_total)
  stage2_test = summary(stage2_test)
  
  stage2X1 = list(theta = stage2_test$coefficients[2,1],
                  vartheta = stage2_test$coefficients[2,2]^2,
                  sigma2 = stage2_test$sigma^2)
  
  FINAL_RESULT = list(gene=gene_exp$gene_name,
                      gene_ind = gene_ind,
                      chr=chr,
                      stage2X1 = stage2X1,
                      fX1 = lm_stage1_X1$fstatistic,
                      lm_stage1_X1 = lm_stage1_X1,
                      n1 = nrow(SNP_BED),
                      n2 = nrow(ukb_pheno)
  )
  real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
}

save(real_data_result,
     file = paste("~/TWAS_testing/results/power/gtex/TWAS_X1_chr",
                  chr,".Rdata",sep=""))
