get_nearest_gene <- function(gene_exp, gene_ind){
  start_pos = gene_exp$V2[gene_ind]
  end_pos = gene_exp$V3[gene_ind]
  start_to_end = abs(gene_exp$V2 - end_pos)
  end_to_start = abs(gene_exp$V3 - start_pos)
  dist = pmin(start_to_end, end_to_start)
  temp = sort(dist, index.return = T)
  if(temp$ix[1]==gene_ind){
    return(temp$ix[2])
  }else{return(temp$ix[1])}
}

get_snp <- function(gene_exp, gene_ind, chr, IGAP_summary, Index_Mat){
  gene_name = as.character(gene_exp$Symbol[gene_ind])
  start = floor(gene_exp[gene_ind,3]/1000) - 100
  end = ceiling(gene_exp[gene_ind,4]/1000) + 100
  
  plink_command = paste0("../twas/plink --bfile ~/TWAS_testing/data/ADNI_SNP/GATK_chr",chr,
                         " --chr ",chr," --from-kb ",
                         start," --to-kb ",end,
                         " --geno 0 --maf 0.05 --hwe 0.001  ",
                         " --make-bed --out ~/TWAS_testing/data/For_Individual_Genes2/",
                         gene_name)
  system(plink_command)
  
  snp = tryCatch(read_plink(paste("~/TWAS_testing/data/For_Individual_Genes2/",gene_name,sep="")),
                 error=function(e){cat("ERROR :",
                                       conditionMessage(e),
                                       "\n")})
  if(is.null(snp))
  {
    return(NULL)
  }
  
  remove_command = paste("rm ~/TWAS_testing/data/For_Individual_Genes2/",
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
    return(NULL)
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
    return(NULL)
  }
  ind = which(is.element(colnames(SNP_BED),colnames(cor_bed)))
  SNP_BED = SNP_BED[,ind]
  snp_bim_igap = snp_bim_igap[ind,]
  
  return(list(SNP_BED = SNP_BED, snp_bim_igap = snp_bim_igap))
}

reg_out_cov <- function(gene_exp, gene_ind, cov5table, SNP, snp_bim){
  # cov5table: covariates
  X = scale(as.numeric(gene_exp[gene_ind,-(1:4)]))
  X1 = X

  lm_RegCov = lm(X1 ~ cov5table$AGE + cov5table$PTGENDER + 
                   cov5table$PTEDUCAT + cov5table$ICV + 
                   cov5table$PTHAND)
  X1 = lm_RegCov$residuals
  X1 = scale(X1)
  if(ncol(SNP) > 50)
  {
    ind = order(abs(cor(X1,SNP)),decreasing = T)[1:50]
    ind = sort(ind)
    SNP = SNP[,ind]
    snp_bim = snp_bim[ind,]
  }
  SNP = scale(SNP)
  return(list(X1 = X1, SNP_BED = SNP, snp_bim_igap = snp_bim))
}

remove_na_effects <- function(X1, SNP, IGAP_cor){
  lm_stage1_X1 = lm(X1 ~ SNP)
  hatbetaX1 = lm_stage1_X1$coefficients[-1]
  
  na_ind = which(!is.na(hatbetaX1))
  SNP = SNP[,na_ind]
  IGAP_cor = IGAP_cor[na_ind,]
  
  return(list(SNP_BED = SNP,
              IGAP_r = IGAP_cor))
}

fit_stage1 <- function(X1, SNP){
  m = 
    step(lm(X1~.,data = data.frame(X1,SNP)),direction = "backward", trace = -1)
  lm_stage1_X1 = summary(m)
  hatbetaX1 = rep(0,ncol(SNP))
  coef_X1 = lm_stage1_X1$coefficients
  name_coef_X1 = substr(rownames(coef_X1),2,30)
  name_SNP_BED = colnames(SNP)
  for(beta_ind in 1:nrow(coef_X1))
  {
    ii = which(name_SNP_BED == name_coef_X1[beta_ind])
    hatbetaX1[ii] = coef_X1[beta_ind,1]
  }
  
  
  if(sum(abs(hatbetaX1)) == 0 )
  {
    return(NULL)
  }
  
  return(list(hatbetaX1 = hatbetaX1, lm_stage1_X1 = m))
}
