setwd("~/TWAS_testing/Rcode")
source("./power/multivariate/power_functions.R")
library(dplyr)
library(tidyr)
library(data.table)

gene_gtf = rtracklayer::import('/home/panwei/shared/GTEx-v8/gencode.v26.GRCh38.genes.gtf')
gene_gtf = as.data.frame(gene_gtf) %>% select(gene_id,gene_name) %>% distinct()
GENE_EXP = fread("~/deepRIV/UKB/data/GTEx/Whole_Blood.v8.normalized_expression.bed.gz")
colnames(GENE_EXP)[1] = 'chrs'
GENE_EXP$chrs = as.numeric(gsub("[^0-9.-]", "", GENE_EXP$chrs))
idx = match(GENE_EXP$gene_id, gene_gtf$gene_id)
GENE_EXP$gene_name = gene_gtf$gene_name[idx]
GENE_EXP = as.data.frame(na.omit(GENE_EXP))

gene_name = GENE_EXP[,c('chrs','gene_name')]
colnames(gene_name) = c('chr', 'Symbol')
gene_name$Symbol = as.character(gene_name$Symbol)
all_result = NULL
all_gene = NULL
all_sig = NULL
s2_factor = 2
r2_cutoff = 0.05

for(chr in 1:22)
{
  load(paste0("~/TWAS_testing/results/power/gtex/TWAS_X1_chr",chr,".Rdata"))
  for(i in 1:length(real_data_result))
  {
    if(is.null(real_data_result[i]$FINAL_RESULT))
    {
      all_result = rbind(all_result,
                         rep(100,11))
      next()
    }
    fX1 = real_data_result[i]$FINAL_RESULT$fX1
    
    if(is.null(fX1)){
      all_result = rbind(all_result,
                         rep(100,11))
      next()
    }
    
    s1_sample = real_data_result[i]$FINAL_RESULT$n1
    s2_sample = round(real_data_result[i]$FINAL_RESULT$n2 * s2_factor)
    print(c(real_data_result[i]$FINAL_RESULT$n1, s2_sample))
    
    stage1_X1 = pf(fX1[1],fX1[2],fX1[3],lower.tail = F)
    sigmaX1 = real_data_result[i]$FINAL_RESULT$stage2X1$sigma2

    # t statistics
    tX1 = real_data_result[i]$FINAL_RESULT$stage2X1$theta / sqrt(real_data_result[i]$FINAL_RESULT$stage2X1$vartheta)
    stage2X1 = 2*pt(abs(tX1), s2_sample-2, lower.tail = F)
    
    X1_theta = real_data_result[i]$FINAL_RESULT$stage2X1$theta
    X1_vartheta = real_data_result[i]$FINAL_RESULT$stage2X1$vartheta
    X1_stage1_err_var = sum(real_data_result[i]$FINAL_RESULT$lm_stage1_X1$residuals^2)/(669)
    X1_stage2_err_var = real_data_result[i]$FINAL_RESULT$stage2X1$sigma2
    X1_r2 = real_data_result[i]$FINAL_RESULT$lm_stage1_X1$r.squared
    
    powerX1 = 
      power_fn(s1_sample,s2_sample,
               X1_stage1_err_var, 
               X1_stage2_err_var, 
               X1_theta,
               X1_vartheta*(real_data_result[i]$FINAL_RESULT$n2 - 1)/(s2_sample - 1))
    X1_vartheta_fml = powerX1$vartheta
    
    all_result = rbind(all_result,
                       c(chr,stage1_X1, stage2X1,
                         X1_r2, X1_theta, X1_vartheta_fml, X1_vartheta, X1_stage1_err_var, X1_stage2_err_var,
                         real_data_result[i]$FINAL_RESULT$n1, real_data_result[i]$FINAL_RESULT$n2))
  }
  gene_chr = which(gene_name$chr == chr)
  gene_chr = gene_name[gene_chr,]
  all_gene = rbind(all_gene,gene_chr)
}
colnames(all_result) = c("chr","stage1_X1", "X1_p",
                         "X1_r2", "X1_theta", "X1_vartheta_fml", "X1_vartheta", "X1_stage1_err_var", "X1_stage2_err_var",
                         "n1", "n2")
all_result = cbind(all_gene,
                   all_result[,-1])

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


all_sig = cbind(all_sig, X1_power_fml, X1_power)

write.table(all_sig,
            paste0("~/TWAS_testing/results/power_curve/gtex/n2/power_table_n2_",s2_factor,"_all.txt"),
            row.names = F)


