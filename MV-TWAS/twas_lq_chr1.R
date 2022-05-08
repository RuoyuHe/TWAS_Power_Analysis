setwd("~/TWAS_testing/Rcode/power/")
library(plink2R)
library(data.table)
source("./allele_qc.R")
source("./stage2.R")
source("./power_functions.R")

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
  
  gene_name = as.character(gene_exp$Symbol[gene_ind])
  start = floor(gene_exp[gene_ind,3]/1000) - 100
  end = ceiling(gene_exp[gene_ind,4]/1000) + 100
  
  plink_command = paste0("../twas/plink --bfile ~/TWAS_testing/data/ADNI_SNP/GATK_chr",chr,
                        " --chr ",chr," --from-kb ",
                        start," --to-kb ",end,
                        " --geno 0 --maf 0.05 --hwe 0.001  ",
                        " --make-bed --out ~/TWAS_testing/data/For_Individual_Genes/",
                        gene_name)
  system(plink_command)
  
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
  
  ### overlap between SNP and IGAP
  colnames(snp$bed) = snp$bim$V4
  
  colnames(snp$bim)[4:6] = c("Position","A1","A2")
  snp_bim_igap = merge(snp$bim,IGAP_summary,by="Position")
  remove_flip = allele.qc(snp_bim_igap$A1,snp_bim_igap$A2,
                          snp_bim_igap$Effect_allele,
                          snp_bim_igap$Non_Effect_allele)
  snp_bim_igap$Beta[which(remove_flip$flip)] = 
    -snp_bim_igap$Beta[which(remove_flip$flip)]
  snp_bim_igap$JansenCor[which(remove_flip$flip)] = 
    -snp_bim_igap$JansenCor[which(remove_flip$flip)]
  
  snp_bim_igap = snp_bim_igap[remove_flip$keep,]
  
  snp_bed_ind = NULL
  snp_bed_colname = colnames(snp$bed)
  for(i in 1:nrow(snp_bim_igap))
  {
    snp_bed_ind = c(snp_bed_ind,which(snp_bed_colname == snp_bim_igap$Position[i]))
  }
  SNP_BED = snp$bed[Index_Mat$fam_ind,snp_bed_ind]
  
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
  SNP_BED = SNP_BED[,ind]
  snp_bim_igap = snp_bim_igap[ind,]
  
  ### regress gene exp on covariates
  
  X = scale(as.numeric(gene_exp[gene_ind,-(1:4)]))
  X1 = X
  X2 = X^2
  X2_ori = X2
  
  lm_RegCov_X1 = lm(X1 ~ cov5table$AGE + cov5table$PTGENDER + 
                   cov5table$PTEDUCAT + cov5table$ICV + 
                   cov5table$PTHAND)
  X1 = lm_RegCov_X1$residuals
  
  lm_RegCov_X2 = lm(X2 ~ cov5table$AGE + cov5table$PTGENDER + 
                   cov5table$PTEDUCAT + cov5table$ICV + 
                   cov5table$PTHAND)
  X2 = lm_RegCov_X2$residuals
  
  X1 = scale(X1)
  X2 = scale(X2)
  
  if(ncol(SNP_BED) > 50)
  {
    ind1 = order(abs(cor(X1,SNP_BED)),decreasing = T)[1:50]
    ind2 = order(abs(cor(X2,SNP_BED)),decreasing = T)[1:50]
    ind = union(ind1,ind2)
    
    ind = sort(ind)
    SNP_BED = SNP_BED[,ind]
    snp_bim_igap = snp_bim_igap[ind,]
  }
  
  SNP_BED = scale(SNP_BED)
  
  ### correlation between SNP and Y
  IGAP_r = snp_bim_igap$Beta / sqrt(54162-2+snp_bim_igap$Beta^2)
  #IGAP_r = snp_bim_igap$Beta / sqrt((54162-2)*snp_bim_igap$SE^2+snp_bim_igap$Beta^2)
  Jansen_r = snp_bim_igap$JansenCor
  ###
  
  lm_stage1_X1 = lm(X1 ~ SNP_BED)
  lm_stage1_X2 = lm(X2 ~ SNP_BED)
  
  hatbetaX1 = lm_stage1_X1$coefficients[-1]
  hatbetaX2 = lm_stage1_X2$coefficients[-1]
  
  na_ind = which(!is.na(hatbetaX1))
  SNP_BED = SNP_BED[,na_ind]
  IGAP_r = IGAP_r[na_ind]
  Jansen_r = Jansen_r[na_ind]
  
  lm_stage1_X1 = 
    step(lm(X1~.,data = data.frame(X1,SNP_BED)),direction = "backward", trace = -1)
  lm_stage1_X1 = summary(lm_stage1_X1)
  hatbetaX1 = rep(0,ncol(SNP_BED))
  coef_X1 = lm_stage1_X1$coefficients
  name_coef_X1 = substr(rownames(coef_X1),2,30)
  name_SNP_BED = colnames(SNP_BED)
  for(beta_ind in 1:nrow(coef_X1))
  {
    ii = which(name_SNP_BED == name_coef_X1[beta_ind])
    hatbetaX1[ii] = coef_X1[beta_ind,1]
  }
  lm_stage1_X1$hatbeta = hatbetaX1
  
  lm_stage1_X2 = 
    step(lm(X2~.,data = data.frame(X2,SNP_BED)),direction = "backward", trace = -1)
  lm_stage1_X2 = summary(lm_stage1_X2)
  hatbetaX2 = rep(0,ncol(SNP_BED))
  coef_X2 = lm_stage1_X2$coefficients
  name_coef_X2 = substr(rownames(coef_X2),2,30)
  name_SNP_BED = colnames(SNP_BED)
  for(beta_ind in 1:nrow(coef_X2))
  {
    ii = which(name_SNP_BED == name_coef_X2[beta_ind])
    hatbetaX2[ii] = coef_X2[beta_ind,1]
  }
  lm_stage1_X2$hatbeta = hatbetaX2
  
  if((sum(abs(hatbetaX1)) + sum(abs(hatbetaX2))) == 0 )
  {
    FINAL_RESULT = NULL
    real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
    
    next()
  }
  
  ### Model 4.2 (see notes)
  W = data.frame(cbind(AGE = cov5table$AGE, GENDER = cov5table$PTGENDER, 
                   EDU = cov5table$PTEDUCAT, ICV =  cov5table$ICV, 
                   HAND = cov5table$PTHAND))
  W$AGE = as.numeric(W$AGE)
  W$EDU = as.numeric(W$EDU)
  W$ICV = as.numeric(W$ICV)
  m.4.2_X1 = step(lm(X~., data = data.frame(X,W,SNP_BED)), direction = "backward",
              scope = list(lower = formula(lm(X~., data = data.frame(X,W)))))
  m.4.2_X2 = step(lm(X2_ori~., data = data.frame(X2_ori,W, SNP_BED)), direction = "backward",
                  scope = list(lower = formula(lm(X2_ori~., data = data.frame(X2_ori,W)))))
  if(length(m.4.2_X1$coefficients)>6){
    smy_m.4.2_X1 = summary(m.4.2_X1)
    m.4.2_X1_design = model.matrix(m.4.2_X1)
    m.4.2_X1_sigma = sum(smy_m.4.2_X1$residuals^2)/smy_m.4.2_X1$df[2]
    m.4.2_squared_X1 = 
      m.4.2_X1_design[,1:6]%*%m.4.2_X1$coefficients[1:6] + 
      (as.matrix(m.4.2_X1_design[,-(1:6)])%*%as.matrix(m.4.2_X1$coefficients[-(1:6)]))^2 + 
      m.4.2_X1_sigma
    
    m.4.2_squared_X1_rsq = var(m.4.2_squared_X1)/var(X2_ori)
    df = smy_m.4.2_X1$df[1]-6
    if (choose(df,2)+df+6<711){
      m.4.2_squared_X1_adj_rsq = 1 - (1 - m.4.2_squared_X1_rsq)*710/(711-choose(df,2)-df-6)
    } else{m.4.2_squared_X1_adj_rsq = NULL}
    
  } else{
    m.4.2_squared_X1_rsq = NULL
    m.4.2_squared_X1_adj_rsq = NULL
  }
  
  rsq.4.2 = list(rsq = summary(m.4.2_X2)$r.squared,
                 adj_rsq = summary(m.4.2_X2)$adj.r.squared,
                 rsq_squared_X1 = m.4.2_squared_X1_rsq[1],
                 adj_rsq_squared_X1 = m.4.2_squared_X1_adj_rsq[1])
  
  
  ### Model 4.3 (see notes on Notability)
  X2_4.3 = X1^2
  m.4.3_X2 = step(lm(X2_4.3~.,data = data.frame(X2_4.3,SNP_BED)),direction = "backward", trace = -1)
  ### get the coef of m.4.3_X2 for stage 2
  m.4.3_X2 = summary(m.4.3_X2)
  m.4.3_hatbetaX2 = rep(0,ncol(SNP_BED))
  m.4.3_coef_X2 = m.4.3_X2$coefficients
  m.4.3_name_coef_X2 = substr(rownames(m.4.3_coef_X2),2,30)
  name_SNP_BED = colnames(SNP_BED)
  for(beta_ind in 1:nrow(m.4.3_coef_X2))
  {
    ii = which(name_SNP_BED == m.4.3_name_coef_X2[beta_ind])
    m.4.3_hatbetaX2[ii] = m.4.3_coef_X2[beta_ind,1]
  }
  m.4.3_X2$hatbeta = m.4.3_hatbetaX2
  
  m.4.3_squared_X1 = SNP_BED%*%hatbetaX1 + sum(lm_stage1_X1$residuals^2)/lm_stage1_X1$df[2]
  m.4.3_squared_X1_rsq = var(m.4.3_squared_X1)/var(X2_4.3)
  df = sum(hatbetaX1!=0)
  if(choose(df,2)+df<711){
    m.4.3_squared_X1_adj_rsq = 1-(1-m.4.3_squared_X1_rsq)*710/(710-choose(df,2)-df)
  } else{m.4.3_squared_X1_adj_rsq=NULL}
  rsq.4.3 = list(rsq = m.4.3_X2$r.squared,
                 adj_rsq = m.4.3_X2$adj.r.squared,
                 rsq_squared_X1 = m.4.3_squared_X1_rsq[1],
                 adj_rsq_squared_X1 = m.4.3_squared_X1_adj_rsq[1])
  
  stage1_squared_X1 = list(rsq.4.2 = rsq.4.2,
                           rsq.4.3 = rsq.4.3)
  
  ### IGAP
  stage2X1 = 
    stage2(betahat = hatbetaX1,
           G = SNP_BED,
           corY = t(t(IGAP_r)),
           n = 54162)
  
  stage2X2 = 
    stage2(betahat = hatbetaX2,
           G = SNP_BED,
           corY = t(t(IGAP_r)),
           n = 54162)
  
  stage2both = 
    stage2(betahat = cbind(hatbetaX1,hatbetaX2),
           G = SNP_BED,
           corY = t(t(IGAP_r)),
           n = 54162)
  
  if(sum(m.4.3_hatbetaX2!=0)>0){
    m.4.3_stage2X2 = 
      stage2(betahat = m.4.3_hatbetaX2,
             G = SNP_BED,
             corY = t(t(IGAP_r)),
             n = 54162)
    m.4.3_stage2both = 
      stage2(betahat = cbind(hatbetaX1,m.4.3_hatbetaX2),
             G = SNP_BED,
             corY = t(t(IGAP_r)),
             n = 54162)
  } else{
    m.4.3_stage2X2 = NULL
    m.4.3_stage2both = NULL
    }
  
  FINAL_RESULT = list(stage2X1 = stage2X1,
                      stage2X2 = stage2X2,
                      stage2both = stage2both,
                      m.4.3_stage2X2 = m.4.3_stage2X2,
                      m.4.3_stage2both =  m.4.3_stage2both,
                      fX1 = lm_stage1_X1$fstatistic,
                      fX2 = lm_stage1_X2$fstatistic,
                      lm_stage1_X1 = lm_stage1_X1,
                      lm_stage1_X2 = lm_stage1_X2,
                      m.4.3_X2 = m.4.3_X2,
                      stage1_squared_X1 = stage1_squared_X1,
                      chr = chr,
                      gene_name = gene_exp$Symbol[gene_ind]
                      
                      
  )
  
  real_data_result = c(real_data_result,list(FINAL_RESULT = FINAL_RESULT))
  
}
save(real_data_result,
     file = paste("~/TWAS_testing/results/power/multivariate/X1_X2/TWAS_X1_chr",
                  chr,".Rdata",sep=""))
