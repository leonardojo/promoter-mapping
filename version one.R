library(tidyverse)
options(stringsAsFactors = F)

table <- read.table('test1.bed.txt', header = F, sep = "\t")
gene.list <- unique(table$V4)
names(table)[4] <- "Genes"
head(table)



ggplot(data = subset(table, Genes == "GENE A")) + 
  geom_point(data = subset(table, V12 == 0.1), aes(x=V8, y= 0.2, fill = V10), shape = 25, size = 3) + 
  geom_point(data = subset(table, V12 == -0.1), aes(x=V9, y= -0.2, fill = V10), shape = 24, size = 3) + 
  xlim(as.numeric(table[2,3])-500, as.numeric(table[2,3])) + 
  ylim(-0.5,0.5) +
  theme_bw() +
  coord_fixed(ratio = 90) +
  facet_wrap("Genes", ncol=1) +
  geom_segment(aes(x=table[1,2],xend=table[1,3],y=0,yend=0))



for (i in gene.list){
ggplot(data = subset(table, Genes == i)) + 
  geom_point(data = subset(table, V12 == 0.1), aes(x=V8, y=V12, fill = V10), shape = 25, size = 2) + 
  geom_point(data = subset(table, V12 == -0.1), aes(x=V9, y=V12, fill = V10), shape = 24, size = 2) + 
  xlim(as.numeric(table[2,3])-500, as.numeric(table[2,3])) + 
  ylim(-2,2) +
  theme_bw() +
  coord_fixed(ratio = 0.2) +
}
  
ggsave(paste(i,".jpeg", sep=""))           
  