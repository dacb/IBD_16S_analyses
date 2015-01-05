library(vegan)

md <- read.delim("sample_info.germ_free.xls", header=T, row.names=1, sep='\t');
head(md)

d <- read.delim("12.persistent_otus.xls", header=T, row.names=1, sep="\t")

# extract just the fecal transplant samples
cols <- row.names(md)
cols
raw_otus <- d[, cols]

# prepare the matrix for processing
sample_otus <- t(raw_otus)
head(sample_otus)


sample_full <- data.frame(sample_otus, cage = md$cage, mouse = paste(md$cage, md$animal, sep=""), week = md$week, gender = md$gender, diet = md$diet)

# drop no rnalater
sample_full[row.names(sample_full) == 'D_donor_nRL', ] <- NA
sample_full[row.names(sample_full) == 'C_donor_nRL', ] <- NA

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

pdf("16.figures.gf.pdf");
plot(mds$points, type = "n")
#, ylim=c(-.5, .5), xlim=c(-.5,.5))
points(mds$points[sample_full$diet == 'control', ], cex = 3.0, pch=21, col="black", bg="white")
points(mds$points[sample_full$diet == 'HVD', ], cex = 3.0, pch=21, col="black", bg="brown")
points(mds$points[sample_full$gender == 'female', ], cex = 2.0, pch=21, col="black", bg="magenta")
points(mds$points[sample_full$gender == 'male', ], cex = 2.0, pch=21, col="black", bg="cyan")
#points(mds$points[sample_full$cage == 1, ], cex = 1.2, pch=21, col="black", bg="green")
#points(mds$points[sample_full$cage == 2, ], cex = 1.2, pch=21, col="black", bg="brown")
#points(mds$points[sample_full$cage == 3, ], cex = 1.2, pch=21, col="black", bg="orange")
points(mds$points[sample_full$week == 1, ], cex = 1.0, pch=21, col="black", bg="white")
points(mds$points[sample_full$week == 2, ], cex = 1.0, pch=21, col="black", bg="grey")
points(mds$points[sample_full$week == 4, ], cex = 1.0, pch=21, col="black", bg="black")
#points(x, y, cex = 2.0, pch=21, col="black", bg="black")
#points(mds$points[sample_full, cex = 2.0, pch=21, col="black", bg="black")
#text(mds$points[sample_full$gender == 'female', ] - 0.06, labels = rownames(sample_full)[sample_full$gender == 'female'], col="magenta", cex=.6)
#text(mds$points[sample_full$gender == 'male', ] - 0.06, labels = rownames(sample_full)[sample_full$gender == 'male'], col="cyan", cex=.6)

#plot(ord_o, type="n")
#text(ord_o, display = "species", cex=.6)
#points(ord_o$CCA$wa[sample_full$diet == 'control diet', ], cex = 0.8, pch=21, col="black", bg="magenta")
#points(ord_o$CCA$wa[sample_full$diet == 'HVD diet', ], cex = 0.8, pch=21, col="black", bg="cyan")
#text(ord_o$CCA$wa[sample_full$diet == 'control diet', ] - 0.06, labels = rownames(sample_full)[sample_full$diet == 'control diet'], col="magenta", cex=.6)
#text(ord_o$CCA$wa[sample_full$diet == 'HVD diet', ] - 0.06, labels = rownames(sample_full)[sample_full$diet == 'HVD diet'], col="cyan", cex=.6)

