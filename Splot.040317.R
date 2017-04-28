library(ggplot2)
options(stringsAsFactors = F)



#####File reading necessary for Splot

#####Make a list for all the BED files for each subregion. 
files <- list.files(pattern = "18-26nt.combined.dblsorted.bed") 

#bed file for stemloop aligned to genome. Used for finding out the size for stemloop(for Splot x axis)
sizeDat <- read.table("All.miRNA.Stemloop.forMirValidation.030716.sorted.bed",sep = "\t", header = F)
#The file for mature miRNA sequence aligned to stemloops to mark miRNA with black triangle. BED file format
miRlist <- read.table("Part1.3PmiRNAmatureSeqAlignToStemloop.042816.txt", sep = "\t", header = F) 
#The file for star strand sequence  aligned to stemloops to make star sequences with black squares.
miRstarList <- read.table("Part1.3PmiRNAStarSeqAlignToStemloop.051316.txt", sep = "\t", header = F)
#Files for the library size (18-26nt) for each subregion. I used to calculate CPM for each sRNAs
srLength <- read.table("Modified.FormiRVal.18-26nt.CombinedSizes.042816.txt", sep = "\t", header = F)


#############################################
#####Generating Splot
################################################

#####List of the stemloop of interest for Splot. 
#### File structure is 1) miRNA mature name, 2) miRNA stemloop name, 3) subregion of interest
###I will send this file to you for your reference. 
listDat <- read.table("Combined.Part1.ForSplot.CPMtableForValidatedMIRwithStemloopAndLibName.052516.txt", sep = "\t", header = T)

#make a new directly for Splot
dir.create("PlotOut.Part1.090916", showWarnings = F)

#Here, I am selecting the BED files that I use for Splots based ont the listDat
unique_file <- unique(listDat[,3])

for(i in c(1:length(unique_file))){
  ###Select stemloop that are in the subregion 
  dat_list <- subset(listDat,listDat[,3]==unique_file[i])
  ###Getting the BED file from the subregion. See "ExampleBEDfile.txt" file for the BED file example
  selected_file <- paste(dat_list[1,3],".qual.clipped.noRDNAorTDNA.Gmax_189.Aligned.18-26nt.combined.dblsorted.bed", sep = "")
  dat <- read.table(selected_file, sep = "\t", header = F)
  MIRname <- dat_list[,2] ####get a list of stemloop 
  samples <- dat_list[1,3] ###get a subergion name
  libSize <- srLength[srLength$V1==samples,2] ####Extract the library size for the subregion 
  mirName <- dat_list[,1] ###get a list of miRNA mature name
  
  for(j in c(1:length(MIRname))){
    name <- MIRname[j] ###Select a miRNA stemloop
    name.miR= mirName[j] ####Select a miRNA mature name
    stemloop <- subset(dat,dat[,1]==name) ###Subset the bed file based on a stemloop name
    MIRsize <- sizeDat[sizeDat$V4==name,3]-sizeDat[sizeDat$V4==name,2] ####finding the length for stemloop
    MIRsize <- MIRsize[1]
    mir<- mirName[j] ###get a miRNA name from the miRNA name list
    
    ###Making a summary table that gives the read counts for each sRNAs.  
    ut <-data.frame(table(stemloop[,c(2,3,6)])) 
    ###Eliminated the sRNA reads with 0 counts
    ut <- subset(ut,ut[,4]>0)
    
    ###I have to add the following scripts to transform tha data.frame vectors in the table to numeric or character
    ut_2 <- as.numeric(as.character(ut$V2))
    ut[,1] <- ut_2
    ut_3 <- as.numeric(as.character(ut$V3))
    ut[,2] <- ut_3
    ut_4 <- as.character(ut$V6)
    ut[,3] <- ut_4
    
    ###Calculate CPM for each read counts
    ut_freq<- ut[,4]/libSize*1000000
    ut[,4] <-  ut_freq
    ####Get the sRNA size
    ut["size"] <-ut$V3-ut$V2
    
    #get negative number if the reads are located on the negative strand
    for(k in c(1:length(ut[,4]))){
      if(ut[k,3]=="-"){
        ut[k,4]=-abs(ut[k,4])
      }
    }
    
    ###Get max and minimum counts  
    maxCount <- max(ut[,4])
    minCount <- min(ut[,4])
    ###Sort the table based on the read counts to get TOP1 and TOP2
    sorted <- sort(ut[,4], decreasing = T)
    top2count <- sorted[2]
    
    #Find out miRNA position in stemloop
    miRBeg <- miRlist[miRlist$V1== name & miRlist$V4==name.miR,2]
    miREnd <- miRlist[miRlist$V1 == name & miRlist$V4==name.miR,3]
    
    #Find out miRNA star position in stemloop
    miRstarBeg <- miRstarList[miRstarList$V1 ==name & miRstarList$V4==name.miR,2]
    miRstarEnd <- miRstarList[miRstarList$V1 ==name & miRstarList$V4==name.miR,3]
    
    miRpos <- ut[ut$V2==miRBeg&ut$V3==miREnd,]
    miRstarpos <- ut[ut$V2==miRstarBeg&ut$V3==miRstarEnd,]
    top1pos <- ut[ut$Freq==maxCount,]
    top2pos <- ut[ut$Freq==top2count,]
    if(is.na(top2pos$Freq)){
      top2pos=top1pos
    }    
    #plot to the designated directory
    filename <- paste(name,"_in_",samples,".pdf", sep = "")
    out <- paste("PlotOut.Part1.090916", filename, sep="/")
    pdf(file=out, width = 5,height = 5)
    print(ggplot(ut,aes(x=ut$V2,y=ut$Freq))+geom_point(
      aes(color=factor(size)))+xlim(0,MIRsize)+ylim(
        c(-abs(minCount)-1,abs(maxCount)+1))+geom_point(
          data=miRpos,aes(x=miRpos$V2,y=miRpos$Freq),
          shape = 17, size = 3)+geom_point(
            data=miRstarpos,aes(x=miRstarpos$V2,y=miRstarpos$Freq),
            shape = 15, size = 3)+labs(title = paste(mir,"-",samples,sep=""),x="Stemloop Position", y="CPM")+geom_text(
              data=top1pos, aes(x=top1pos$V2,y=top1pos$Freq, label = "top1"))+geom_text(
                data=top2pos, aes(x=top2pos$V2,y=top2pos$Freq, label = "top2")))
    dev.off()
  }
}
