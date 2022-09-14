setwd("~/Dropbox/UMN/Research_Projects/GWAS_testing/TWAS/Rcode/power/power_curve/multivariate")
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(dplyr)
library(tidyr)
library(data.table)


read_data <- function(){
  x1x2 <<- fread("./IGAP_X1X2/chisq_power_table_all.txt",header = T, data.table = F)
  x1near <<- fread("./IGAP_nearest/chisq_power_table_all.txt",header = T, data.table = F)
  uni <<- fread("../SampleSize/IGAP_only/power_table_1_all.txt",header = T, data.table = F)
  
  common_symbols = intersect(x1x2$Symbol, x1near$Symbol)
  common_symbols = intersect(common_symbols, uni$Symbol)
  x1x2 <<- x1x2 %>% filter(x1x2$Symbol %in% common_symbols)
  x1near <<- x1near %>% filter(x1near$Symbol %in% common_symbols)
  uni <<- uni %>% filter(uni$Symbol %in% common_symbols)
  print(all(x1near$Symbol==uni$Symbol))
}

cutoff = 0.01
cbPalette <- c("#D55E00", "#56B4E9", "#009E73", 
               "#0072B2", "#CC79A7") 

#------------------------ compare X1
read_data()
idx = which((uni$X1_p<=cutoff) & (x1x2$X1_p<=cutoff) & (x1near$X1_p<=cutoff))
a_x1x2 = x1x2[idx,]
a_x1near = x1near[idx,]
a_uni = uni[idx,]

# UV-L vs LQ-L vs MV-Target, p-value <= 0.01
a_x = 1:length(idx)
a_s = sort(a_uni$X1_power_fml, index.return = T)
p1 = ggplot() + 
  xlab("genes") +
  ylab("power") +
  ggtitle("A") +
  geom_line(aes(x = a_x, y = a_uni$X1_power_fml[a_s$ix], linetype = "UV-L", colour = "UV-L")) + 
  geom_line(aes(x = a_x, y = a_x1x2$X1_power_fml[a_s$ix], linetype = " TWAS-L",  colour = " TWAS-L")) +
  geom_line(aes(x = a_x, y = a_x1near$X1_power_fml[a_s$ix], linetype = "MV-Target", colour = "MV-Target")) +
  scale_linetype(name = '') + 
  scale_colour_manual(values = cbPalette, name = '') + 
  theme(legend.position = c(0.10, 0.75),
        legend.text = element_text(size=12),
        legend.background = element_rect(fill=alpha('white',0)))

tmp = x1x2$X1_power_fml - uni$X1_power_fml
idx = sort(tmp, decreasing = T, index.return = T)$ix
x1x2[idx[1:5],]
tmp[idx[1:5]]


tmp2 = x1near$X1_power_fml - uni$X1_power_fml
idx = sort(tmp2, decreasing = T, index.return = T)$ix
x1near[idx[1:5],]
tmp2[idx[1:5]]

#------------------------- compare X2
read_data()
idx = which((uni$X2_p<=cutoff) & (x1x2$X2_p<=cutoff) & (x1near$X2_p<=cutoff))
b_x1x2 = x1x2[idx,]
b_x1near = x1near[idx,]
b_uni = uni[idx,]

#UV-Q vs LQ-Q vs MV-Alternative, p-value <= 0.01
b_x = 1:length(idx)
b_s = sort(b_uni$X2_power_fml, index.return = T)
p2 = ggplot() + 
  xlab("genes") +
  ylab("power") +
  ggtitle("B") +
  geom_line(aes(x = b_x, y = b_uni$X2_power_fml[b_s$ix], linetype = "UV-Q", colour = "UV-Q")) + 
  geom_line(aes(x = b_x, y = b_x1x2$X2_power_fml[b_s$ix], linetype = " TWAS-Q",  colour = " TWAS-Q")) +
  geom_line(aes(x = b_x, y = b_x1near$X2_power_fml[b_s$ix], linetype = "MV-Alternative",  colour = "MV-Alternative")) +
  scale_linetype(name = '') + 
  scale_colour_manual(values = cbPalette, name = '') + 
  theme(legend.position = c(0.11, 0.75),
        legend.text = element_text(size=12),
        legend.background = element_rect(fill=alpha('white',0)))

tmp = x1x2$X2_power_fml - uni$X2_power_fml
idx = sort(tmp, decreasing = T, index.return = T)$ix
x1x2[idx[1:5],]
tmp[idx[1:5]]


tmp2 = x1near$X2_power_fml - uni$X2_power_fml
idx = sort(tmp2, decreasing = T, index.return = T)$ix
x1near[idx[1:5],]
tmp2[idx[1:5]]

#------------------------- comapre X1 to both
read_data()
idx = which((uni$X1_p<=cutoff) & (x1x2$stage2_p<=cutoff) & (x1near$stage2_p<=cutoff))
c_x1x2 = x1x2[idx,]
c_x1near = x1near[idx,]
c_uni = uni[idx,]

# UV-L vs LQ-LQ vs MV-Joint, p-value <= 0.01

arrows <- 
  tibble(
    x1 = c(49, 55 + 3, 6, 41, 32),
    x2 = c(49, 53.2, 6, 44.7, 30),
    y1 = c(0.965, 0.23, 0.225, 0.84, 0.18), 
    y2 = c(0.85, 0.034, 0.12, 0.99, 0.11)
  )

c_x = 1:length(idx)
c_s = sort(c_uni$X1_power_fml, index.return = T)
p3 = ggplot() + 
  xlab("genes") +
  ylab("power") +
  ggtitle("C") +
  geom_line(aes(x = c_x, y = c_uni$X1_power_fml[c_s$ix], linetype = "UV-L", colour = "UV-L")) + 
  geom_line(aes(x = c_x, y = c_x1x2$chisq_both_power_fml[c_s$ix], linetype = " TWAS-LQ", colour = " TWAS-LQ")) +
  geom_line(aes(x = c_x, y = c_x1near$chisq_both_power_fml[c_s$ix], linetype = "MV-Joint", colour = "MV-Joint")) +
  scale_linetype(name = '') + 
  scale_colour_manual(values = cbPalette, name = '') + 
  annotate(geom="text", x=49, y=1, label="HLA-DRB5",
           color="red") +
  annotate(geom="text", x=53 + 5, y=0.25, label="GLT8D2",
           color="red") +
  annotate(geom="text", x=6, y=0.27, label="RNASE3",
           color="red") +
  annotate(geom="text", x=41, y=0.8, label="HLA-DRB1",
           color="red") +
  annotate(geom="text", x=32, y=0.22, label="RNASE2",
           color="red") +
  geom_curve(
    data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2),
    arrow = arrow(length = unit(0.08, "inch")), size = 0.5,
    color = "gray20", curvature = 0) +
  theme(legend.position = c(0.10, 0.75),
        legend.text = element_text(size=12),
        legend.background = element_rect(fill=alpha('white',0)))

tmp = x1x2$chisq_both_power_fml - uni$X1_power_fml
idx = sort(tmp, decreasing = T, index.return = T)$ix
x1x2[idx[1:5],]
tmp[idx[1:5]]
uni$X1_power_fml[idx[1:5]]

tmp2 = x1near$chisq_both_power_fml - uni$X1_power_fml
idx = sort(tmp2, decreasing = T, index.return = T)$ix
x1near[idx[1:5],]
tmp2[idx[1:5]]

---------------------# compare X2 to both
read_data()
idx = which((uni$X2_p<=cutoff) & (x1x2$stage2_p<=cutoff) & (x1near$stage2_p<=cutoff))
d_x1x2 = x1x2[idx,]
d_x1near = x1near[idx,]
d_uni = uni[idx,]

#UV-Q vs LQ-LQ vs MV-Joint, p-value <= 0.01
d_x = 1:length(idx)
d_s = sort(d_uni$X2_power_fml, index.return = T)
p4 = ggplot() + 
  xlab("genes") +
  ylab("power") +
  ggtitle("D") +
  geom_line(aes(x = d_x, y = d_uni$X2_power_fm[d_s$ix], linetype = "UV-Q", colour = "UV-Q")) + 
  geom_line(aes(x = d_x, y = d_x1x2$chisq_both_power_fml[d_s$ix], linetype = " TWAS-LQ",  colour = " TWAS-LQ")) +
  geom_line(aes(x = d_x, y = d_x1near$chisq_both_power_fml[d_s$ix], linetype = "MV-Joint",  colour = "MV-Joint")) +
  scale_linetype(name = '') + 
  scale_colour_manual(values = cbPalette, name = '') + 
  theme(legend.position = c(0.10, 0.75),
        legend.text = element_text(size=12),
        legend.background = element_rect(fill=alpha('white',0)))

tmp = x1x2$chisq_both_power_fml - uni$X2_power_fml
idx = sort(tmp, decreasing = T, index.return = T)$ix
x1x2[idx[1:5],]
tmp[idx[1:5]]
  

tmp2 = x1near$chisq_both_power_fml - uni$X2_power_fml
idx = sort(tmp2, decreasing = T, index.return = T)$ix
x1near[idx[1:5],]
tmp2[idx[1:5]]

ggarrange(p1, p2, p3, p4, ncol=1, nrow=4, common.legend = F)

