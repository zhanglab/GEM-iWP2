wd <- commandArgs(trailingOnly = TRUE)
setwd(wd)

library(ggplot2)
library(ggalt)
library(patchwork)

df <- read.csv('./random-effic-values.tsv', header=F, sep='\t')
colnames(df) <- c('temp', 'metric', 'val')
df$temp <- factor(df$temp, levels = c('04C', '15C', '20C'), labels = c('4°C', '15°C', '20°C'), ordered = T)
df$metric <- factor(df$metric, levels = c('CUE', 'ATPP', 'ATPC'), labels = c('Carbon Use Efficiency (CUE)', 'ATP Produced / GlcNac', 'ATP Consumed / gDW Biomass'), ordered = T)

g <- ggplot(df, aes(temp, val,fill=temp))
g+geom_boxplot() + facet_wrap(~metric, scales='free')+
  scale_fill_manual(values=c("blue","gold","firebrick1")) +
  theme_bw() + xlab('Temperature') + ylab('Value') +
  theme(strip.text.x = element_text(face='bold', size=12),
        axis.title.x = element_text(face='bold', size=12),
        axis.title.y = element_text(face='bold', size=12),
        axis.text = element_text(face='bold', size=10))

ggsave('./Figure2.pdf', width=10, height=4)
