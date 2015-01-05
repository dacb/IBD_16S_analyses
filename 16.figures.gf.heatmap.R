#source("ggplot-heatmap.R")
source("heatmap.3.R")
#library(gplots)
library(vegan)
library(RColorBrewer)
#library(gplots)
#library(DESeq2)
#library(fdrtool)
#library(MASS)
library(ggplot2)
#library(reshape)

md <- read.delim("sample_info.germ_free.xls", header=T, row.names=1, sep='\t');
head(md)

d <- read.table("12.persistent_otus.xls", header=T, row.names=1, sep="\t")

# extract just the germ free fecal transplant samples
cols <-row.names(md)
cols[cols == "XDonor"] <- "Donor"
cols
raw_otus <- d[, cols]
head(raw_otus)

# remove rows with all 0s
raw_otus <- raw_otus[apply(raw_otus[,-1], 1, function(x) !all(x==0)),]

sample_read_sums <- read.table("13.summarize_persistent_otus.gf.read_sum.xls", header=T, sep="\t")
sample_read_sums
dim(sample_read_sums)

# normalize the raw_otu matrix to % of reads
otus_pcnt <- scale(raw_otus, center=FALSE, scale=t(sample_read_sums))
head(otus_pcnt)

write.table(otus_pcnt, file="16.figures.otus_percent.tab", quote=F, sep="\t")

# normalize the otu maxtrix by the maximum number of reads
max_seqs <- max(sample_read_sums)
otus <- scale(raw_otus, center=FALSE, scale=t(sample_read_sums))
otus <- otus * max_seqs
colSums(otus)
head(otus)

log_otus <- log(otus)
# artificial baseline of the 0 entries
log_otus[log_otus == -Inf] <- NA
head(log_otus)

#ind <- apply(log_otus, 1, function(x) all(is.na(x)))
#log_otus <- log_otus[ !ind, ]

write.table(log_otus, file="16.figures.log_otus.tab", quote=F, sep="\t")

pdf("./16.figures.gf.heatmap.pdf")

heatmap.3(as.matrix(log_otus), as.matrix(otus_pcnt), key=T, distfun=vegdist, col=brewer.pal(9, "YlGnBu"), hclustfun = function(x) hclust(x,method = 'average'), scale='none', na.rm=T, na.color="white", density.info='none', trace='none', KeyValueName="ln Reads", symkey=F, symbreaks=F, sepcolor="lightgrey", rowsep=c(1:100), colsep=c(1:100))

dev.off()


