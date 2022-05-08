setwd("~/Dropbox/UMN/Research_Projects/GWAS_testing/TWAS/Rcode/power/power_curve/gtex/stage2_results_for_r2")
source("~/Dropbox/UMN/Research_Projects/GWAS_testing/TWAS/Rcode/power/multivariate/power_functions_multivariate.R")
library(dplyr)
library(tidyr)
library(data.table)
library(RSQLite)
library(stringr)

herit <- dbConnect(RSQLite::SQLite(), 
                   "~/Dropbox/UMN/Research_Projects/GWAS_testing/TWAS/Rcode/power/multivariate/gtex/GenArchDB-master/genarch.db")
sqlcommand = paste0('SELECT gene, tissue, h2, h2_ci FROM results WHERE tissue = "WholeBlood_TW";')
res <- dbSendQuery(herit, sqlcommand)
herit = dbFetch(res)
rm(res)
tmp_lower = numeric(nrow(herit))
tmp_upper = numeric(nrow(herit))
a = strsplit(herit$h2_ci, "-")
for(i in 1:length(a)){
  if(length(a[[i]])==3){
    tmp_lower[i] = 0
    tmp_upper[i] = as.numeric(a[[i]][3])
  }else if(length(a[[i]])==2){
    tmp_lower[i] = as.numeric(a[[i]][1])
    tmp_upper[i] = as.numeric(a[[i]][2])
  }
}
herit$lower = tmp_lower
herit$upper = tmp_upper

all_result = NULL
all_gene = NULL
all_sig = NULL
gene_name = NULL

for(chr in 1:22)
{
  load(paste0("TWAS_X1_chr",chr,".Rdata"))
  for(i in 1:length(real_data_result))
  {
    if(is.null(real_data_result[i]$FINAL_RESULT)) next()

    fX1 = real_data_result[i]$FINAL_RESULT$fX1
    if(is.null(fX1)) next()
    
    current_gene = real_data_result[i]$FINAL_RESULT$gene
    if(!(current_gene %in% herit$gene)) next()
    
    X1_r2 = real_data_result[i]$FINAL_RESULT$lm_stage1_X1$r.squared
    current_herit = herit[which(herit$gene==current_gene)[1],c('h2','lower','upper')]
    if(current_herit$h2 < X1_r2) next()
    
    X1_theta = real_data_result[i]$FINAL_RESULT$stage2X1$theta
    X1_vartheta = real_data_result[i]$FINAL_RESULT$stage2X1$vartheta
    X1_stage1_err_var = sum(real_data_result[i]$FINAL_RESULT$lm_stage1_X1$residuals^2)/(669)
    X1_stage2_err_var = real_data_result[i]$FINAL_RESULT$stage2X1$sigma2
    s1_sample = real_data_result[i]$FINAL_RESULT$n1
    s2_sample = real_data_result[i]$FINAL_RESULT$n2
    
    h2_r2_diff = current_herit$h2 - X1_r2
    h2_X1_stage1_err_var = X1_stage1_err_var - h2_r2_diff
    
    upper_r2_diff = current_herit$upper - X1_r2
    upper_X1_stage1_err_var = X1_stage1_err_var - upper_r2_diff
    
    if(current_herit$lower > 0){
      lower_r2_diff = current_herit$lower - X1_r2
      lower_X1_stage1_err_var = X1_stage1_err_var - lower_r2_diff
    }else{
      lower_X1_stage1_err_var = X1_stage1_err_var
    }
    
    powerX1 = 
      power_fn(s1_sample,s2_sample,
               X1_stage1_err_var, 
               X1_stage2_err_var, 
               X1_theta,
               X1_vartheta)
    X1_vartheta_fml = powerX1$vartheta
    
    h2_powerX1 = 
      power_fn(s1_sample,s2_sample,
               h2_X1_stage1_err_var, 
               X1_stage2_err_var, 
               X1_theta,
               X1_vartheta)
    h2_X1_vartheta_fml = h2_powerX1$vartheta
    
    upper_powerX1 = 
      power_fn(s1_sample,s2_sample,
               upper_X1_stage1_err_var, 
               X1_stage2_err_var, 
               X1_theta,
               X1_vartheta)
    upper_X1_vartheta_fml = upper_powerX1$vartheta
    
    lower_powerX1 = 
      power_fn(s1_sample,s2_sample,
               lower_X1_stage1_err_var, 
               X1_stage2_err_var, 
               X1_theta,
               X1_vartheta)
    lower_X1_vartheta_fml = lower_powerX1$vartheta
    
    stage1_X1 = pf(fX1[1],fX1[2],fX1[3],lower.tail = F)
    sigmaX1 = real_data_result[i]$FINAL_RESULT$stage2X1$sigma2

    # t statistics
    tX1 = real_data_result[i]$FINAL_RESULT$stage2X1$theta / sqrt(real_data_result[i]$FINAL_RESULT$stage2X1$vartheta)
    stage2X1 = 2*pt(abs(tX1), s2_sample-2, lower.tail = F)
    
    gene_name = c(gene_name, current_gene)
    all_result = rbind(all_result,
                       c(chr, stage1_X1, stage2X1,
                         X1_r2, current_herit$h2, current_herit$upper, current_herit$lower,
                         X1_theta, X1_vartheta_fml, 
                         h2_X1_vartheta_fml, upper_X1_vartheta_fml, lower_X1_vartheta_fml,
                         X1_vartheta, X1_stage1_err_var, X1_stage2_err_var,
                         real_data_result[i]$FINAL_RESULT$n1, real_data_result[i]$FINAL_RESULT$n2))
  }
}
all_result = as.data.frame(all_result)
colnames(all_result) = c("chr", "stage1_X1", "X1_p",
                         "X1_r2", "h2", "h2_upper", "h2_lower",
                         "X1_theta", "X1_vartheta_fml", 
                         "h2_X1_vartheta_fml", "upper_X1_vartheta_fml", "lower_X1_vartheta_fml",
                         "X1_vartheta", "X1_stage1_err_var", "X1_stage2_err_var",
                         "n1", "n2")
all_result$gene_name = gene_name


valid_genes = sum((all_result$stage1_X1<0.001))
print(paste("valid genes:",valid_genes))

all_result = all_result[all_result$stage1_X1<0.001,]

# stage2all = which((all_result$X1_p < 0.01))

# all_sig = all_result[stage2all,]
all_sig = all_result

X1_power_fml = mapply(get_power, 
                      all_sig$X1_theta, all_sig$X1_vartheta_fml, rep(0.05, dim(all_sig)[1]),
                      rep(valid_genes, dim(all_sig)[1]))
X1_power = mapply(get_power, 
                  all_sig$X1_theta, all_sig$X1_vartheta, rep(0.05, dim(all_sig)[1]),
                  rep(valid_genes, dim(all_sig)[1]))
h2_X1_power_fml = mapply(get_power, 
                      all_sig$X1_theta, all_sig$h2_X1_vartheta_fml, rep(0.05, dim(all_sig)[1]),
                      rep(valid_genes, dim(all_sig)[1]))
upper_X1_power_fml = mapply(get_power, 
                         all_sig$X1_theta, all_sig$upper_X1_vartheta_fml, rep(0.05, dim(all_sig)[1]),
                         rep(valid_genes, dim(all_sig)[1]))
lower_X1_power_fml = mapply(get_power, 
                            all_sig$X1_theta, all_sig$lower_X1_vartheta_fml, rep(0.05, dim(all_sig)[1]),
                            rep(valid_genes, dim(all_sig)[1]))

all_sig = cbind(all_sig, X1_power_fml, X1_power,
                h2_X1_power_fml, upper_X1_power_fml, lower_X1_power_fml)

write.table(all_sig,
            paste0("~/Dropbox/UMN/Research_Projects/GWAS_testing/TWAS/Rcode/power/power_curve/gtex/combined_power_table_r2.txt"),
            row.names = F)


