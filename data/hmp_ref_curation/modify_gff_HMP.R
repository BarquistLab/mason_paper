library(dplyr)

setwd("~/Documents/mason_paper/data/hmp_ref_curation")

args <- c("./start_regs_new.gff", "./genelengths.tsv", )


# import gff and fasta lengths:
gff <-  read.delim(args[1], header = FALSE) 

fastas <- read.delim(args[2], header = FALSE, row.names = 1)


# fet length of fasta for each gff entry:
gff$RL <- fastas[gff$V1,]

errorgenes <- gff[gff$RL < gff$V5 ,]$V3

#change gff if RL is bigger than end:
gff$V5 <- ifelse(gff$RL < gff$V5, gff$RL, gff$V5)

# delete rows which are equal:
gff <- gff[!(duplicated(gff$V4) & duplicated(gff$V5) & duplicated(gff$V3)),]
gff <- gff[!(gff$V5<gff$V4),]

# write gff that is modified:
write.table(gff[!(gff$RL < gff$V5 | gff$V4 < 1),-10],
          file = "start_regs_new_mod.gff",  
          sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


