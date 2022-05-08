library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(latex2exp)

s1_factor = c(0,1,2,3)
a = fread(paste0("~/TWAS_testing/results/power_curve/r2/p_001/power_table_",0,".txt"),
          data.table = F)
for (i in s1_factor){
  b = fread(paste0("~/TWAS_testing/results/power_curve/r2/p_001/power_table_",
                   i,".txt"),data.table=F)
  print(all(a$Symbol==b$Symbol))
  a = cbind(a,b$X1_power_fml)
  colnames(a)[ncol(a)] = paste0("X1_power_fml_",i)
}
fwrite(a,"~/TWAS_testing/results/power_curve/r2/p_001/combined_p001.txt")


########### Sample Size #############
##### n1
setwd("~/Dropbox/UMN/Research_Projects/GWAS_testing/TWAS/Rcode/power/power_curve/gtex")
a = fread("./combined_power_table_n1.txt",header = T, data.table = F)
a1 = a %>% filter(a$X1_p <= 0.05/nrow(a))
a = a %>% filter(a$X1_p <= 0.01)
x = 1:nrow(a)
s = sort(a$X1_power, index.return = T)

cc <- scales::seq_gradient_pal("#63fa05", "#057ffa", "Lab")(seq(0,1,length.out=5))
cbPalette <- c("#D55E00", "#009E73", "#56B4E9",
               "#0072B2", "#6c40bd", "#CC79A7")
p1 = ggplot() + 
  xlab("genes") +
  ylab("power") +
  ggtitle("A") +
  geom_line(aes(x = x, y = a$X1_power_fml[s$ix], linetype = ' 1 x', colour = ' 1 x')) + 
  geom_line(aes(x = x, y = a$X1_power_fml_2[s$ix], linetype = ' 2 x', colour = ' 2 x')) +
  geom_line(aes(x = x, y = a$X1_power_fml_5[s$ix], linetype = ' 5 x', colour = ' 5 x')) +
  geom_line(aes(x = x, y = a$X1_power_fml_10[s$ix], linetype = '10 x', colour = '10 x')) +
  geom_line(aes(x = x, y = a$X1_power_fml_50[s$ix], linetype = '50 x', colour = '50 x')) +
  geom_line(aes(x = x, y = a$X1_power[s$ix], linetype = 'Inf', colour = 'Inf')) +
  geom_vline(xintercept=nrow(a) - nrow(a1), linetype="dashed", 
             color = "red", size=0.5) + 
  annotate(geom="text", x=nrow(a) - nrow(a1) + 3, y=0.125, label="B",
             color="red") +
  scale_linetype(name = '') + 
  scale_colour_manual(values = cbPalette, name = '') + 
  theme(legend.text = element_text(size=15),
        legend.position = c(0.05, 0.87))

  
if (F){scale_linetype_manual(name = '', values = c('x1','x2','x5','x10','Inf'),
                             labels = c("x1", "x2" , " x5" , 
                                        "x10" , unname(TeX("$\\infty$"))))}

# x1 = (nrow(a)-nrow(a1)+1):nrow(a)
x1 = 1:nrow(a1)
s1 = sort(a1$X1_power, index.return = T)
p2 = ggplot() + 
  xlab("genes") +
  ylab("power") +
  ggtitle("A") +
  geom_line(aes(x = x1, y = a1$X1_power_fml[s1$ix], linetype = ' 1 x', colour = ' 1 x')) + 
  geom_line(aes(x = x1, y = a1$X1_power_fml_2[s1$ix], linetype = ' 2 x', colour = ' 2 x')) +
  geom_line(aes(x = x1, y = a1$X1_power_fml_5[s1$ix], linetype = ' 5 x', colour = ' 5 x')) +
  geom_line(aes(x = x1, y = a1$X1_power_fml_10[s1$ix], linetype = '10 x', colour = '10 x')) +
  geom_line(aes(x = x1, y = a1$X1_power_fml_50[s1$ix], linetype = '50 x', colour = '50 x')) +
  geom_line(aes(x = x1, y = a1$X1_power[s1$ix], linetype = 'Inf', colour = 'Inf')) +
  scale_linetype(name = '') + 
  scale_colour_manual(values = cbPalette, name = '') + 
  theme(legend.text = element_text(size=15),
        legend.position = c(0.05, 0.87))

ggarrange(p1, p2, ncol=1, nrow=2, common.legend = TRUE, legend="right")

##### n2
setwd("~/Dropbox/UMN/Research_Projects/GWAS_testing/TWAS/Rcode/power/power_curve/gtex")
a = fread("./combined_power_table_n2.txt",header = T, data.table = F)

a1 = a %>% filter(a$X1_p <= 0.05/nrow(a))
a = a %>% filter(a$X1_p <= 0.01)
x = 1:nrow(a)
s = sort(a$X1_power, index.return = T)

cc <- scales::seq_gradient_pal("#63fa05", "#057ffa", "Lab")(seq(0,1,length.out=5))
cbPalette <- c("#D55E00", "#009E73", "#56B4E9",
               "#0072B2", "#6c40bd", "#CC79A7")
p1 = ggplot() + 
  xlab("genes") +
  ylab("power") +
  ggtitle("A") +
  geom_line(aes(x = x, y = a$X1_power_fml[s$ix], linetype = '1 x', colour = '1 x')) + 
  geom_line(aes(x = x, y = a$X1_power_fml_12[s$ix], linetype = '1.2 x', colour = '1.2 x')) +
  geom_line(aes(x = x, y = a$X1_power_fml_15[s$ix], linetype = '1.5 x', colour = '1.5 x')) +
  geom_line(aes(x = x, y = a$X1_power_fml_2[s$ix], linetype = '2 x', colour = '2 x')) +
  geom_vline(xintercept=nrow(a) - nrow(a1), linetype="dashed", 
             color = "red", size=0.5) + 
  annotate(geom="text", x=nrow(a) - nrow(a1) + 3, y=0.125, label="B",
           color="red") +
  scale_linetype(name = '') + 
  scale_colour_manual(values = cbPalette, name = '') + 
  theme(legend.text = element_text(size=15),
        legend.position = c(0.05, 0.87))


# x1 = (nrow(a)-nrow(a1)+1):nrow(a)
x1 = 1:nrow(a1)
s1 = sort(a1$X1_power, index.return = T)
p2 = ggplot() + 
  xlab("genes") +
  ylab("power") +
  ggtitle("A") +
  geom_line(aes(x = x1, y = a1$X1_power_fml[s1$ix], linetype = '1 x', colour = '1 x')) + 
  geom_line(aes(x = x1, y = a1$X1_power_fml_12[s1$ix], linetype = '1.2 x', colour = '1.2 x')) +
  geom_line(aes(x = x1, y = a1$X1_power_fml_15[s1$ix], linetype = '1.5 x', colour = '1.5 x')) +
  geom_line(aes(x = x1, y = a1$X1_power_fml_2[s1$ix], linetype = '2 x', colour = '2 x')) +
  scale_linetype(name = '') + 
  scale_colour_manual(values = cbPalette, name = '') + 
  theme(legend.text = element_text(size=15),
        legend.position = c(0.05, 0.87))

ggarrange(p1, p2, ncol=1, nrow=2, common.legend = TRUE, legend="right")



########### r2 #############
setwd("~/Dropbox/UMN/Research_Projects/GWAS_testing/TWAS/Rcode/power/power_curve/gtex")
library(grid)
library(gtable)
library(ggpubr)
library(latex2exp)
a = fread("combined_power_table_r2.txt",header = T, data.table = F)

a = a %>% filter(X1_p <= 0.05/nrow(a))

cbPalette <- c("#D55E00", "#56B4E9") 

x = 1:nrow(a)
s = sort(a$upper_X1_power_fml, index.return = T)
a = a[s$ix,]
a$x = x
p1 = ggplot(data=a, aes(x=x, y=h2_X1_power_fml, colour = 'heritability')) + 
  geom_point() +
  geom_line() + 
  geom_ribbon(aes(ymin=lower_X1_power_fml, ymax=upper_X1_power_fml), linetype=2, alpha=0.1) + 
  geom_line(aes(x = x, y = X1_power_fml, colour = "baseline")) +
  xlab("genes") +
  ylab("power") +
  ggtitle("A") +
  scale_colour_manual(values = cbPalette, name = '') +
  theme(legend.text = element_text(size=15),
        legend.position = c(0.1, 0.75),
        legend.background = element_rect(fill=alpha('white',0)))


p2 = ggplot(data=a, aes(x=x, y=h2, colour = 'heritability')) + 
  xlab("genes") +
  ylab(unname(TeX("$R^2$"))) + 
  geom_point() +
  geom_line() + 
  geom_ribbon(aes(ymin=h2_lower, ymax=h2_upper), linetype=2, alpha=0.1) + 
  geom_line(aes(x = x, y = X1_r2, colour = "baseline")) +
  ggtitle("B") +
  scale_linetype(name = '') +
  scale_colour_manual(values = cbPalette, name = '') + 
  theme(legend.text = element_text(size=15),
        legend.position = c(0.1, 0.75),
        legend.background = element_rect(fill=alpha('white',0)))
  
p3 = ggplot(a, aes(x=x, y=abs(X1_theta[s$ix]))) +
  xlab("ordered genes") +
  ylab(unname(TeX("$|\\hat{\\theta}_i|$"))) + 
  ggtitle("C") +
  geom_point(size=1, colour = "black", alpha = 0.4)

p4 = ggplot(a0, aes(x=x, y=X1_vartheta_fml[s$ix])) +
  xlab("ordered genes") +
  ylab(unname(TeX("$\\widehat{Var}_a(\\hat{\\theta}_i)$"))) + 
  ggtitle("D") +
  geom_point(size=1, colour = "black", alpha = 0.4, shape = 18)

p5 = ggplot(data = NULL, aes(x=x, y = a3$X1_vartheta_fml[s$ix]/a0$X1_vartheta_fml[s$ix])) +
  xlab("ordered genes") +
  ylab(unname(TeX("ratio of $\\widehat{Var}_a(\\hat{\\theta}_i)$"))) + 
  ggtitle("E") +
  geom_point(size=1, colour = "black", alpha = 0.4, shape = 23)

ggarrange(p1, p2, p3, p4, p5, ncol=1, nrow=5, common.legend = TRUE, legend="right")
ggarrange(p1, p2, ncol=1, nrow=2, common.legend = T)

g1 <- ggplotGrob(p1)
g1 <- gtable_add_cols(g1, unit(0,"mm")) # add a column for missing legend
g2 <- ggplotGrob(p2)
g3 <- ggplotGrob(p3)
g1$widths <- unit.pmax(g1$widths, g2$widths, g3$widths)
g2$widths <- unit.pmax(g1$widths, g2$widths, g3$widths)
g3$widths <- unit.pmax(g1$widths, g2$widths, g3$widths)

# stack them afterwards
g <- rbind(g1, g2, g3, size="first") # stack the two plots

g$layout[grepl("guide", g$layout$name),c("t","b")] <- c(1,nrow(g))
grid.newpage()
grid.draw(g)
########### rf #############
setwd("~/Dropbox/UMN/Research_Projects/GWAS_testing/TWAS/Rcode/power/power_curve/rf")
a = fread("power_table_p01.txt",header = T, data.table = F)
x = 1:nrow(a)
s = sort(a$X1_power_rf, index.return = T)
ggplot() + 
  xlab("ordered genes") +
  ylab("power") +
  ggtitle("p-value <= 0.01") +
  geom_line(aes(x = x, y = a$X1_power_fml[s$ix], colour = "  baseline")) + 
  geom_line(aes(x = x, y = a$X1_power_rf[s$ix], colour = " rf"))



varplot = ggplot(a, aes(x = X1_r2, y = X1_r2_rfvar)) + 
  geom_point(color = "red", shape = 23)+
  xlab("lm r2") +
  ylab("rf r2") +
  ggtitle("p <= 0.1, var(X_hat)/var(X)") +
  geom_abline(slope = 1)

corplot = ggplot(a, aes(x = X1_r2, y = X1_r2_rfcor)) + 
  geom_point(color = "red", shape = 23)+
  xlab("lm r2") +
  ylab("rf r2") +
  ggtitle("p <= 0.1, var(X_hat)/var(X)") +
  geom_abline(slope = 1)
