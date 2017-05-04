# Promoter Mapping v1.0

library(tidyverse)
options(stringsAsFactors = F)

## Organizing the table
table <- read.table('Motifs_Promoters.bed', header = F, sep = "\t")
table <- table[,-7]
names(table)[7:12] <- c("V7","V8","V9","V10",'V11',"V12")
names(table)[4] <- "Genes"
head(table)

## Setting up coordinates for each motif in the promoter based on their position on chromossomes

table$V8 <- table$V8 - table$V2 ## Start
table$V9 <- table$V9 - table$V2 ## Stop
head(table)

## Setting up coordinates for START and STOP for each motif adjusting Strand Orientation, 
## For instance, if a gene is on the negative strand, motifs that are in position 10 will be placed on (PROMOTER SIZE) - 10

table$motif.start[table$V6 == "+"] <- table[table$V6 == "+", 8] ## Gene on positive strand, just copy the start from orignal
table$motif.stop[table$V6 == "+"] <- table[table$V6 == "+", 9] ## Gene on positive strand, just copy the STOP from orignal

table$motif.start[table$V6 == "-"] <- (table[table$V6 == "-", 3]-table[table$V6 == "-", 2]) - table[table$V6 == "-", 9] ## Gene on negative strand, promoter size (V3 - V2) minus stop
table$motif.stop[table$V6 == "-"] <- (table[table$V6 == "-", 3]-table[table$V6 == "-", 2]) - table[table$V6 == "-", 8] ## Gene on negative strand, promoter size (V3 - V2) minus start

## Setting up the middle point of the motif
table$motif.middle <- table$motif.start + as.numeric((table$motif.stop-table$motif.start)/2)

head(table)

## Subset table from a gene list
gene.list <- read.table("EM LEARBZAB DT.txt", header = F)
table.s <- merge(table, gene.list, by.x = "Genes", by.y = "V1")

## Preparing for ploting

##Setting up scale points
promoter.size <- as.numeric(table[1,3] - table[1,2]-1) ## They all need to be the same
point <- data.frame(seq.int(0,promoter.size, by = 100))
names(point) <- "X"

##Setting up image size
h <- n_distinct(table.s$Genes)/4

ggplot(data = table.s) + 
  
  ## Motifs in the positive strand of genes in positive strand
  geom_point(data = subset(table.s, (V6 == "+") & (V12 == "+") & (V10 != "TSS")), aes(x=motif.middle, y= 0.25, fill = V10), shape = 25, size = 2) + 
  ## Motifs in the negative strand of genes in positive strand
  geom_point(data = subset(table.s, (V6 == "+") & (V12 == "-") & (V10 != "TSS")), aes(x=motif.middle, y= -0.25, fill = V10), shape = 24, size = 2, show.legend = FALSE) + 
  ## Motifs in the positive strand of genes in negative strand
  geom_point(data = subset(table.s, (V6 == "-") & (V12 == "+") & (V10 != "TSS")), aes(x=motif.middle, y= -0.25, fill = V10), shape = 24, size = 2, show.legend = FALSE) +
  ## Motifs in the negative strand of genes in negative strand
  geom_point(data = subset(table.s, (V6 == "-") & (V12 == "-") & (V10 != "TSS")), aes(x=motif.middle, y= 0.25, fill = V10), shape = 25, size = 2) +
  
  ##TSS
  geom_point(data = subset(table.s, V10 == "TSS"), aes(x=motif.middle, y= 0, colour= "TSS"), fill = "grey20", shape = 18, size = 2) + 
  
  geom_point(data = point, aes(x= X, y =0.0), fill= "grey20", shape = 108, size = 2) +
  
  xlim(0,promoter.size) + 
  ylim(-0.80,0.80) +
  
  theme_bw() +
  coord_fixed(ratio = 100) +
  
  facet_wrap("Genes", ncol= 3) +
  geom_segment(aes(x=0,xend=promoter.size,y=0,yend=0)) +
  scale_fill_manual(values=c("red2","lightgreen","dodgerblue")) +
  scale_colour_manual(values=c("grey20")) +
  theme(legend.position = "top", legend.title = element_blank(),
        axis.title.y=element_blank(), axis.title.x = element_blank(), axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  
  ggsave(file = "LABA.jpeg", dpi = 600, width = 5, height = as.numeric(h), units = "in", limitsize = F)


dev.off()


