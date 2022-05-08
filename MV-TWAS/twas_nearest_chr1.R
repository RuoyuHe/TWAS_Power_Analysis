setwd("~/TWAS_testing/Rcode/power/")
library(plink2R)
library(data.table)
source("./allele_qc.R")
source("./stage2.R")
source("./multivariate/X1_nearest_gene/stage1.R")

chr = 1

### index matrix
Index_Mat = read.table("~/TWAS_testing/data/GenerateIndex_Feb4.txt", header = T)

### IGAP
IGAP_summary = fread(paste0("~/TWAS_testing/data/22traits/Imputed_IGAP_by_Chr/all_IGAP_chr", chr),
                          header = T, data.table = F)
R2_filter = which(IGAP_summary$Var>=0.3)
IGAP_summary = IGAP_summary[R2_filter,]
IGAP_summary = cbind(chr,IGAP_summary[,c(2,1,3,4,5)],1,1)
colnames(IGAP_summary) = c("Chromosome","Position","MarkerName","Effect_allele",    
                             "Non_Effect_allele","Beta","SE","Pvalue"  )
### Jansen GWAS
JansenGWAS = fread(paste0("~/TWAS_testing/data/22traits/Jansen/Jansen_chr",
                              chr,".txt"), header = T, data.table = F)
cat(dim(JansenGWAS))
JansenGWAS = JansenGWAS[JansenGWAS$Nsum > 300000,]
cat(dim(JansenGWAS))

### Merge Two GWAS
JansenCor = JansenGWAS$Z / sqrt(JansenGWAS$Nsum-2+(JansenGWAS$Z)^2)
# JansenCor = JansenGWAS$BETA / sqrt((JansenGWAS$Nsum-2)*JansenGWAS$SE^2 + (JansenGWAS$BETA)^2)
JansenGWAS = cbind(JansenGWAS[,2:5],JansenCor)
colnames(JansenGWAS)[1:2] = c("Chromosome","Position")
IGAP_summary = merge(IGAP_summary,JansenGWAS,by = c("Chromosome","Position"))
IGAP_summary$Effect_allele = as.character(IGAP_summary$Effect_allele)
IGAP_summary$Non_Effect_allele = as.character(IGAP_summary$Non_Effect_allele)
IGAP_summary$A1 = as.character(IGAP_summary$A1)
IGAP_summary$A2 = as.character(IGAP_summary$A2)
TwoGWASqc = allele.qc(IGAP_summary$Effect_allele,IGAP_summary$Non_Effect_allele,
                      IGAP_summary$A1,IGAP_summary$A2)
IGAP_summary$JansenCor[which(TwoGWASqc$flip)] = 
  -IGAP_summary$JansenCor[which(TwoGWASqc$flip)]
IGAP_summary = IGAP_summary[TwoGWASqc$keep,]
IGAP_summary = IGAP_summary[,-c(9:10)]

### Gene Expression
gene_exp = fread("~/TWAS_testing/data/gene_exp_Feb4.txt",header = T,data.table = F)
gene_exp = gene_exp[gene_exp$V1 == chr,
                    c(1:4,Index_Mat$cov5_ind+4)]

### Covariates
cov5table = fread("~/TWAS_testing/data/cov5table_Feb4.txt",header = T,data.table=F)
cov5table = cov5table[Index_Mat$cov5_ind,]

real_data_result = list()
# start analysis ----------------------------------------------------------
cat(nrow(gene_exp),"\n")
for(gene_ind in 1:nrow(gene_exp))
{
  cat(gene_ind,"\n")
  
  snp_X1 = get_snp(gene_exp, gene_ind, chr, IGAP_summary, Index_Mat)
  if (is.null(snp_X1)){
    FINAL_RESULT = NULL
    real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    next()
  } else{
    SNP_BED = snp_X1$SNP_BED
    snp_bim_igap = snp_X1$snp_bim_igap
  }
  
  gene_ind_nearest = get_nearest_gene(gene_exp, gene_ind)
  snp_X1_nearest = get_snp(gene_exp, gene_ind_nearest, chr, IGAP_summary, Index_Mat)
  if (is.null(snp_X1_nearest)){
    FINAL_RESULT = NULL
    real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    next()
  } else{
    SNP_BED_nearest = snp_X1_nearest$SNP_BED
    snp_bim_igap_nearest = snp_X1_nearest$snp_bim_igap
  }
  
  ### regress gene exp on covariates
  snp_X1 = reg_out_cov(gene_exp, gene_ind, cov5table, SNP_BED, snp_bim_igap)
  X1 = snp_X1$X1
  SNP_BED = snp_X1$SNP_BED
  snp_bim_igap = snp_X1$snp_bim_igap
  
  snp_X1_nearest = reg_out_cov(gene_exp, gene_ind_nearest, cov5table, 
                               SNP_BED_nearest, snp_bim_igap_nearest)
  X1_nearest = snp_X1_nearest$X1
  SNP_BED_nearest = snp_X1_nearest$SNP_BED
  snp_bim_igap_nearest = snp_X1_nearest$snp_bim_igap
  
  ### correlation between SNP and Y
  IGAP_r = snp_bim_igap$Beta / sqrt(54162-2+snp_bim_igap$Beta^2)
  IGAP_r = as.data.frame(cbind(IGAP_r = IGAP_r, Position = snp_bim_igap$Position))
  #IGAP_r = snp_bim_igap$Beta / sqrt((54162-2)*snp_bim_igap$SE^2+snp_bim_igap$Beta^2)
  Jansen_r = snp_bim_igap$JansenCor
  Jansen_r = as.data.frame(cbind(Jansen_r = Jansen_r, Position = snp_bim_igap$Position))
  
  IGAP_r_nearest = snp_bim_igap_nearest$Beta / sqrt(54162-2+snp_bim_igap_nearest$Beta^2)
  IGAP_r_nearest = as.data.frame(cbind(IGAP_r = IGAP_r_nearest, 
                                       Position = snp_bim_igap_nearest$Position))
  Jansen_r_nearest = snp_bim_igap_nearest$JansenCor
  Jansen_r_nearest = as.data.frame(cbind(Jansen_r = Jansen_r_nearest,
                                         Position = snp_bim_igap_nearest$Position))
  
  # remove NA effects
  rmNA_X1 = remove_na_effects(X1, SNP_BED, IGAP_r, Jansen_r)
  rmNA_X1_nearest = remove_na_effects(X1_nearest, SNP_BED_nearest, 
                                      IGAP_r_nearest, Jansen_r_nearest)
  
  # backward selection
  stage1X1 = fit_stage1(X1, rmNA_X1$SNP_BED)
  if(is.null(stage1X1)){
    FINAL_RESULT = NULL
    real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    next()
  }
  
  stage1X1_nearest = fit_stage1(X1_nearest, rmNA_X1_nearest$SNP_BED)
  if(is.null(stage1X1_nearest)){
    FINAL_RESULT = NULL
    real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    next()
  }
  
  ####### preparation for stage2both #######
  # combine data, remove duplicate SNPs
  SNP_BED_combined = cbind(rmNA_X1$SNP_BED, rmNA_X1_nearest$SNP_BED)
  idx_dup = which(duplicated(colnames(SNP_BED_combined)))
  if (length(idx_dup)>0){
    SNP_BED_combined = SNP_BED_combined[,-idx_dup]
    
    IGAP_r_combined = rbind(rmNA_X1$IGAP_r, rmNA_X1_nearest$IGAP_r)[-idx_dup,]
    Jansen_r_combined = rbind(rmNA_X1$Jansen_r, rmNA_X1_nearest$Jansen_r)[-idx_dup,]
    temp = as.numeric(colnames(SNP_BED_combined))
  } else{
    IGAP_r_combined = rbind(rmNA_X1$IGAP_r, rmNA_X1_nearest$IGAP_r)
    Jansen_r_combined = rbind(rmNA_X1$Jansen_r, rmNA_X1_nearest$Jansen_r)
    temp = as.numeric(colnames(SNP_BED_combined))
  }
 
  # check all SNPs are aligned
  print(paste0("IGAP_r & SNPs aligned? ", all(IGAP_r_combined$Position==temp)))
  print(paste0("Jansen_r & SNPs aligned? ", all(Jansen_r_combined$Position==temp)))
  
  # align beta with snps
  hatbetaX1 = c(stage1X1$hatbetaX1, 
                rep(0,ncol(SNP_BED_combined) - length(stage1X1$hatbetaX1)))
  hatbetaX1_nearest = rep(0,ncol(SNP_BED_combined))
  idx = match(colnames(rmNA_X1_nearest$SNP_BED),colnames(SNP_BED_combined))
  hatbetaX1_nearest[idx] = stage1X1_nearest$hatbetaX1
  
  # get statistics
  lm_stage1_X1 = summary(stage1X1$lm_stage1_X1)
  lm_stage1_X1_nearest = summary(stage1X1_nearest$lm_stage1_X1)
  
  
  ####### stage2 #######
  ### IGAP
  stage2X1 = 
    stage2(betahat = stage1X1$hatbetaX1,
           G = rmNA_X1$SNP_BED,
           corY = t(t(rmNA_X1$IGAP_r$IGAP_r)),
           n = 54162)
  
  stage2X1_nearest = 
    stage2(betahat = stage1X1_nearest$hatbetaX1,
           G = rmNA_X1_nearest$SNP_BED,
           corY = t(t(rmNA_X1_nearest$IGAP_r$IGAP_r)),
           n = 54162)
  
  stage2both = 
    stage2(betahat = cbind(hatbetaX1,hatbetaX1_nearest),
           G = SNP_BED_combined,
           corY = t(t(IGAP_r_combined$IGAP_r)),
           n = 54162)
  
  ### Jansen
  'stage2X1_Jansen = 
    stage2(betahat = hatbetaX1,
           G = SNP_BED,
           corY = t(t(Jansen_r)),
           n = 455258)
  
  stage2X2_Jansen = 
    stage2(betahat = hatbetaX2,
           G = SNP_BED,
           corY = t(t(Jansen_r)),
           n = 455258)
           
  stage2both_Jansen = 
    stage2(betahat = cbind(hatbetaX1,hatbetaX2),
           G = SNP_BED,
           corY = t(t(Jansen_r)),
           n = 455258)
  
  FINAL_RESULT = list(stage2X1 = stage2X1,
                      stage2X2 = stage2X2,
                      stage2both = stage2both,
                      stage2X1_Jansen = stage2X1_Jansen,
                      stage2X2_Jansen = stage2X2_Jansen,
                      stage2both_Jansen = stage2both_Jansen,
                      fX1 = lm_stage1_X1$fstatistic,
                      fX2 = lm_stage1_X2$fstatistic
  )'
  
  FINAL_RESULT = list(stage2X1 = stage2X1,
                      stage2X1_nearest = stage2X1_nearest,
                      stage2both = stage2both,
                      fX1 = lm_stage1_X1$fstatistic,
                      fX2 = lm_stage1_X1_nearest$fstatistic,
                      lm_stage1_X1 = lm_stage1_X1,
                      lm_stage1_X1_nearest = lm_stage1_X1_nearest
                      
  )
  
  real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
  
}
save(real_data_result,
     file = paste("~/TWAS_testing/results/power/multivariate/X1_nearest_gene/TWAS_X1_chr",
                  chr,".Rdata",sep=""))
