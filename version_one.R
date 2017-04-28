library(tidyverse)
options(stringsAsFactors = F)

## Organizing the table
table <- read.table('test1.bed.txt', header = F, sep = "\t")

table$V8 <- table$V8 - table$V2
table$V9 <- table$V9 - table$V2

names(table)[4] <- "Genes"


table$motif.start[table$V6 == "+"] <- table[table$V6 == "+", 8]
table$motif.stop[table$V6 == "+"] <- table[table$V6 == "+", 9]
table$motif.start[table$V6 == "-"] <- (table[table$V6 == "-", 3]-table[table$V6 == "-", 2]) - table[table$V6 == "-", 9]
table$motif.stop[table$V6 == "-"] <- (table[table$V6 == "-", 3]-table[table$V6 == "-", 2]) - table[table$V6 == "-", 8]

table$motif.middle <- table$motif.start + as.numeric((table$motif.stop-table$motif.start)/2)

head(table)
point <- data.frame(c(1,100,200,300,400,500))
names(point) <- "X"
 
ggplot(data = table) + 
  geom_point(data = subset(table, (V6 == "+") & (V12 == "+")), aes(x=motif.middle, y= 0.15, fill = V10), shape = 25, size = 3) + 
  geom_point(data = subset(table, (V6 == "+") & (V12 == "-")), aes(x=motif.middle, y= -0.15, fill = V10), shape = 24, size = 3, show.legend = FALSE) + 
  geom_point(data = subset(table, (V6 == "-") & (V12 == "+")), aes(x=motif.middle, y= -0.15, fill = V10), shape = 24, size = 3, show.legend = FALSE) +
  geom_point(data = subset(table, (V6 == "-") & (V12 == "-")), aes(x=motif.middle, y= 0.15, fill = V10), shape = 25, size = 3) +
  geom_point(data = point, aes(x= X, y =0.0), fill= "grey20", shape = 108, size = 3) +

  xlim(1,500) + 
  ylim(-0.5,0.5) +
  theme_bw() +
  coord_fixed(ratio = 100) +
  facet_wrap("Genes", ncol= 3) +
  geom_segment(aes(x=1,xend=500,y=0,yend=0)) +
  scale_fill_manual(values=c("red2","lightgreen","dodgerblue", "grey20")) +
  scale_shape_manual(values=c(24,25)) +
  theme(legend.position = "top", legend.title = element_blank(),
        axis.title.y=element_blank(), axis.title.x = element_blank(), axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggsave("Test.jpeg", width = 10, height = 10)


dev.off()
