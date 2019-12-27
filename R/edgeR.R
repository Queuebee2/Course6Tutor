### HAN (c) Todt

library(edgeR)

###############################################################
### read counts
###############################################################

fDir <-  "C:\\Users\\sanne\\OneDrive\\Documenten\\HAN 2\\Blok6\\"
fName <- "WCFS1_cnts.txt"

cnts <- read.delim(paste0(fDir,fName), comment.char="#")
WCFS1_anno <- read.delim2("~/HAN 2/Blok6/WCFS1_anno.txt")
### used for topTags identification
row.names(cnts) <- cnts[,"ID"]

# makes groups names
exp <- c("WCFS1.glc","WCFS1.glc","WCFS1.rib","WCFS1.rib")

# makes a factor of those group names
group <- factor(exp)

# makes a DGE list for easy editing 
# using only the counts( colom 2 to 5) 
# using the factors as groups
y <- DGEList(counts=cnts[,2:5],group=group)

###############################################################
### Normalise counts
### Trimmed mean of M values : remove lowest and highest values
### (percentile) and calculate mean
### 
###############################################################

y <- calcNormFactors(y, method="TMM" )

###############################################################
### Check statistics
###############################################################

print("Count statistics")
print(summary(y$counts))
print(y$samples)

###############################################################
### Create design matrix to let R know what is what
###############################################################

design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
print(design)

###############################################################
### Estimate Dispersion
### method (low-level function) used to estimated the trended dispersions. 
### Possible values are "auto" (default) 
### switch to "bin.spline" method if the number of genes is great than 200 
### and "power" method otherwise)
###############################################################

#y <- estimateDisp(y, design)

y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design, method="power")
y <- estimateGLMTagwiseDisp(y,design)

###############################################################
### Plot results
### MDS :Plot samples on a two-dimensional scatterplot so that distances
### on the plot approximate the typical log2 fold changes between the samples.
### BCV : Plot the genewise biological coefficient of variation (BCV)
### against gene abundance (in log2 counts per million).
###############################################################


plotMDS(y)
plotBCV(y)


###############################################################
### Fit data
### Fit a negative binomial generalized log-linear model
### to the read counts for each gene. Conduct genewise
### statistical tests for a given coefficient or coefficient
### contrast.
###############################################################

fit <- glmFit(y,design)

###############################################################
### Determine fold changes
###############################################################

# Construct the contrast matrix corresponding to specified contrasts of a set of parameters.
mc <- makeContrasts(exp.r=WCFS1.glc-WCFS1.rib, levels=design)
# mc : Contrasts
#Levels      exp.r
#WCFS1.glc     1
#WCFS1.rib    -1

# glmLRT conducts likelihood ratio tests for one or more coefficients in the linear model.
# If contrast is non-null, then the null hypothesis is that the specified contrasts 
# of the coefficients are equal to zero. For example, a contrast of c(0,1,-1), assuming there are three coefficients, would test the hypothesis that the second and third coefficients are equal.

fit <- glmLRT(fit, contrast=mc)

###############################################################
### Print top tags
### Extracts the most differentially expressed genes (or sequence tags) from a test object
### ranked either by p-value or by absolute log-fold-change.
###############################################################

res<-topTags(fit)
print(res)

###############################################################
### maakt fold change plot
###############################################################

plotMD(fit)
summary(decideTests(fit))

###############################################################
### maakt subset van differentialy expressed genes 
### only significant so pvalue lower than 0.05
### no low fold change so lower than -2 or higher than 2
###############################################################

DE_Genes<- fit$table[which(fit$table$PValue < 0.05), names(fit$table) %in% c("logFC","LogCPM","LR","PValue")]
DE_Genes<- DE_Genes[which(DE_Genes$logFC < -2 | DE_Genes$logFC > 2), names(DE_Genes) %in% c("logFC","LR","PValue")]

###############################################################
### joint de annotatie file met de DE_genes
###############################################################

library(dplyr)
names1 <- row.names(DE_Genes)
names2 <- DE_Genes[1]
names <- data.frame(genes = names1, Foldchange = names2 )
comp <- merge(WCFS1_anno,names, by.x = "ORF",by.y = "genes")

##############################################################
### maakt een lijst met alle EC nummers en haalt alles zonder eruit
###############################################################

sub <- select(filter(comp, comp$EC != ""),c("EC"))
comp2 <- semi_join(comp,sub, by = c("EC" = "EC"))

##############################################################
### schrijft data weg naar een file
###############################################################

write.table(sub, file = "sub.txt")
write.csv(comp2, file = "MyData.csv",row.names=FALSE)

