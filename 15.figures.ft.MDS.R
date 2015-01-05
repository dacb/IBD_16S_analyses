library(vegan)

md <- read.delim("sample_info.fecal_trans.xls", header=T, row.names=1, sep='\t');
head(md)

d <- read.delim("12.persistent_otus.xls", header=T, row.names=1, sep="\t")

# extract just the fecal transplant samples
cols <- paste("X", c(row.names(md)), sep="")
cols[cols == "XDonor"] <- "Donor"
cols
raw_otus <- d[, cols]

# now we normalize first, see next section
# prepare the matrix for processing
#sample_otus <- t(raw_otus)
#head(sample_otus)

sample_read_sums <- read.table("13.summarize_persistent_otus.ft.read_sum.xls", header=T, sep="\t")
sample_read_sums
dim(sample_read_sums)

# normalize the otu maxtrix by the maximum number of reads
max_seqs <- max(sample_read_sums)
otus <- scale(raw_otus, center=FALSE, scale=t(sample_read_sums))
otus <- otus * max_seqs
colSums(otus)
head(otus)
sample_otus <- t(otus)
head(sample_otus)

sample_full <- data.frame(sample_otus, cage = md$cage, mouse = paste(md$cage, md$animal, sep=""), week = md$week, gender = md$gender)

mds <- metaMDS(sample_otus, maxit=1000, trymax=100)
mds

#ord_o <- cca(sample_otus ~ cage + mouse + week + gender, data = sample_full)
#ord_o
#summary(ord_o)

#ord <- cca(sample_otus ~ diet + Animal + Week + Gender, data = sample_full)
#anova(ord)
#anova(ord, by="term")
#anova(ord, by="mar")
#anova(ord, by="axis")

#ord_c_c <- cca(sample_otus ~ Condition(diet) + Replicate + Week, data = sample_full)
#anova(ord_c_c)
#anova(ord_c_c, by="term")
#anova(ord_c_c, by="term", strata=diet)
#anova(ord_c_c, by="mar")
#anova(ord_c_c, by="axis")

mds$points
x=mds$points[is.na(sample_full$cage),1]
y=mds$points[is.na(sample_full$cage),2]
mds$points[sample_full$gender == 'female', ]
pdf("15.figures.ft.MDS.pdf");
plot(mds$points, type = "n", xlim=c(-.5, .6), ylim=c(-.5, .6))
points(mds$points[sample_full$gender == 'female', ], cex = 2.0, pch=21, col="black", bg="magenta")
points(mds$points[sample_full$gender == 'male', ], cex = 2.0, pch=21, col="black", bg="cyan")
#points(mds$points[sample_full$cage == 1, ], cex = 1.2, pch=21, col="black", bg="green")
#points(mds$points[sample_full$cage == 2, ], cex = 1.2, pch=21, col="black", bg="brown")
#points(mds$points[sample_full$cage == 3, ], cex = 1.2, pch=21, col="black", bg="orange")
points(mds$points[sample_full$week == 1, ], cex = 1.0, pch=21, col="black", bg="white")
points(mds$points[sample_full$week == 2, ], cex = 1.0, pch=21, col="black", bg="grey")
points(mds$points[sample_full$week == 4, ], cex = 1.0, pch=21, col="black", bg="black")
points(x, y, cex = 2.0, pch=21, col="black", bg="black")
#text(mds$points[sample_full$gender == 'female', ] - 0.06, labels = rownames(sample_full)[sample_full$gender == 'female'], col="magenta", cex=.6)
#text(mds$points[sample_full$gender == 'male', ] - 0.06, labels = rownames(sample_full)[sample_full$gender == 'male'], col="cyan", cex=.6)

#plot(ord_o, type="n")
#text(ord_o, display = "species", cex=.6)
#points(ord_o$CCA$wa[sample_full$diet == 'control diet', ], cex = 0.8, pch=21, col="black", bg="magenta")
#points(ord_o$CCA$wa[sample_full$diet == 'HVD diet', ], cex = 0.8, pch=21, col="black", bg="cyan")
#text(ord_o$CCA$wa[sample_full$diet == 'control diet', ] - 0.06, labels = rownames(sample_full)[sample_full$diet == 'control diet'], col="magenta", cex=.6)
#text(ord_o$CCA$wa[sample_full$diet == 'HVD diet', ] - 0.06, labels = rownames(sample_full)[sample_full$diet == 'HVD diet'], col="cyan", cex=.6)

