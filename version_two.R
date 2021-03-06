# Promoter Mapping v1.3
## install.packages("tidyverse")
library(tidyverse)
library(dplyr)
library(ggplot2)
options(stringsAsFactors = F)

## Organizing the table
## Table Motifs_Promoters.bed are the result of intersect bed 
## of a BED file containing the coordinates of the promoter 600bp (-500,+100)
## and a BED file containing the coordinates of the motifs and TSS
## Files are in ~/BED_Files
## Promoter Chr - Promoter Start - Promoter Stop - Genes - phytozomev10 - Promoter Strand - Motif Chr - Motif Start - Motif Stop - Motif Name - Motid Score - Motif Strand
table <- read.table('Motifs.Promoter.bed', header = F, sep = "\t")
table <- table[,-7] ## Removing column 7, represents the number of each motif

## Adding names to the Promoter Side of the table (P.***)
names(table)[1:6] <- c("P.Chr","P.Start","P.Stop","Genes",'Phyto',"P.Strand")

## Adding names to the Motif Side of the table (M.***)
names(table)[7:12] <- c("M.Chr","M.Start","M.Stop","Motif",'M.Score',"M.Strand")

head(table)

## Setting up coordinates for each motif in the promoter based on their position on chromossomes
table$Coord.Start.temp <- table$M.Start - table$P.Start ## Start
table$Coord.Stop.temp <- table$M.Stop - table$P.Start ## Stop
head(table)

## Setting up coordinates for START and STOP for each motif adjusting Strand Orientation, 
## For instance, if a gene is on the negative strand, motifs that start in position 10 and end in position 15 will be placed on (PROMOTER SIZE) - 15

table$Coord.Start[table$P.Strand == "+"] <- 
  table[table$P.Strand == "+", 13] ## Maintain START if Gene is in the + Strand
table$Coord.Stop[table$P.Strand == "+"] <- 
  table[table$P.Strand == "+", 14] ## Maintain STOP if Gene is in the + Strand

table$Coord.Start[table$P.Strand == "-"] <- 
  (table[table$P.Strand == "-", 3]- table[table$P.Strand == "-", 2]) - table[table$P.Strand == "-", 14] ## Gene on negative strand, promoter size (Column 3 - Column 2) minus Coord.stop
table$Coord.Stop[table$P.Strand == "-"] <- 
  (table[table$P.Strand == "-", 3]- table[table$P.Strand == "-", 2]) - table[table$P.Strand == "-", 13] ## Gene on negative strand, promoter size (V3 - V2) minus start

head(table)

## Removing TEMP Columns
table <- table[-c(13,14)]

## Setting up the middle point of the motif
table$Motif.Middle <- table$Coord.Start + as.numeric((table$Coord.Stop-table$Coord.Start)/2)

head(table)

## Subset table from a gene list
## Gene List is a tab delimited text file with a list of genes that you want to subset
gene.groups <- list.files(pattern = "LABA_DT_LIPID_STORAGE.txt")

for (i in gene.groups){
gene.list <- read.table(paste(i), header = F)

table.subset <- merge(table, gene.list, by.x = "Genes", by.y = "V1")

head(table.subset)

## Preparing for ploting

##Setting up scale points in the map (it will add scale into the maps 100 bp)
promoter.size <- as.numeric(table[1,3] - table[1,2]-1) ## They all need to be the same
point <- data.frame(seq.int(0,promoter.size, by = 100))
names(point) <- "X"

##Do you want to check the plot first? I recommend reducing the number of genes
## table.subset <- table.subset[c(1:30),]

##Setting up image size
h <- n_distinct(table.subset$Genes)/8

## In the future: Order genes based on expression
name <- paste(substr(i,1,nchar(i)-4),"png",sep=".")
a <- ggplot(data = table.subset) + 
  
  ## Motifs in the positive strand of genes in positive strand
  geom_point(data = subset(table.subset, (P.Strand == "+") & (M.Strand == "+") & (Motif != "TSS")), aes(x=Motif.Middle, y= 0.4, fill = Motif), shape = 25, size = 2) + 
  ## Motifs in the negative strand of genes in positive strand
  geom_point(data = subset(table.subset, (P.Strand == "+") & (M.Strand == "-") & (Motif != "TSS")), aes(x=Motif.Middle, y= -0.4, fill = Motif), shape = 24, size = 2, show.legend = FALSE) + 
  ## Motifs in the positive strand of genes in negative strand
  geom_point(data = subset(table.subset, (P.Strand == "-") & (M.Strand == "+") & (Motif != "TSS")), aes(x=Motif.Middle, y= -0.4, fill = Motif), shape = 24, size = 2, show.legend = FALSE) +
  ## Motifs in the negative strand of genes in negative strand
  geom_point(data = subset(table.subset, (P.Strand == "-") & (M.Strand == "-") & (Motif != "TSS")), aes(x=Motif.Middle, y= 0.4, fill = Motif), shape = 25, size = 2, show.legend = FALSE) +
  
  #TSS
  geom_point(data = subset(table.subset, Motif == "TSS"), aes(x=Motif.Middle, y= 0, colour= "TSS"), fill = "grey20", shape = 18, size = 2) + 
  
  # Adding the scale points in the promoter
  geom_point(data = point, aes(x= X, y =0.0), fill= "grey20", shape = 108, size = 2) +
  
  xlim(0,promoter.size) + 
  ylim(-0.80,0.80) +
  
  theme_bw() +
  coord_fixed(ratio = 100) +
  
  facet_wrap("Genes", ncol= 5) +
  geom_segment(aes(x=0,xend=promoter.size,y=0,yend=0)) +
  scale_fill_manual(values=c("green3","red2","greenyellow","dodgerblue")) +
  scale_colour_manual(values=c("grey20")) +
  theme(legend.position = "top", legend.title = element_blank(),
        axis.title.y=element_blank(), axis.title.x = element_blank(), axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

a + 
  ggsave(file = name, dpi = 600, width = 10, height = as.numeric(h)*1.5, units = "in", limitsize = F)

}

