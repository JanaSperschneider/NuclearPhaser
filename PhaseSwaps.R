library(ggplot2)
#---------------------------------
#---------------------------------
#---------------------------------
#---------------------------------
# Change these two line to your directory & contig of interest
#---------------------------------
setwd("R:/NuclearPhaser")
#---------------------------------
contig <- 'tig00000348'
#---------------------------------
#---------------------------------
#---------------------------------
#---------------------------------
df <- read.delim(paste(contig, "_HiC_Contacts.txt", sep=""), sep= "\t", header=FALSE)
head(df)
#---------------------------------
alignments <- read.delim(paste(contig, "_Haplotigs.txt", sep=""), sep= "\t", header=FALSE)
head(alignments)
#---------------------------------
data <- data.frame(x = df$V1,
                   haplotype0 = df$V4, 
                   haplotype1 = df$V5)
head(data)
maximum <- max(df$V2)
#---------------------------------
align_coords <- data.frame(x = alignments$V2,
                           y = alignments$V3,
                           name = alignments$V1,
                           haplotype0 = alignments$V4)
align_coords
#---------------------------------
require(gridExtra)

options(scipen=10000)

p <- ggplot(data) +
  geom_area(aes(x=x, y=haplotype0), fill = "#8c510a", colour = "#8c510a", alpha=0.5) +      
  stat_smooth(geom = 'area', aes(x=x, y=haplotype0), span = 1/3, alpha=0.5, fill = "#01665e") +
  xlab("Position on contig (Mb)") + 
  ylab("% of Hi-C trans reads\n that link to haplotype 0")  +
  scale_x_continuous(breaks = round(seq(0, maximum, by = 500000),0)) + 
  theme_bw(base_size = 48, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle=90, hjust=1)) + 
  geom_segment(data = align_coords, aes(x = x, y = haplotype0, xend = y, yend = haplotype0, colour = "black"), colour = "black",  size = 3,   lineend = "round") +
  ggtitle(contig) +
  scale_x_continuous(breaks = c(seq(from = 0, to = maximum, by = 500000)),
                     labels = c(seq(from = 0, to = maximum, by = 500000))/1000000,
                     expand = c(0, 0), limits = c(0, NA)) 
  

p 
#---------------------------------
png(paste(contig, ".png", sep=""), height = 10, width = 20, units = 'in', res = 300)
p
dev.off()
