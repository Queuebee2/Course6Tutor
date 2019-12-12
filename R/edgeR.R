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

###############################################################
### Create DGEList
###############################################################

exp <- c("WCFS1.glc","WCFS1.glc","WCFS1.rib","WCFS1.rib")

group <- factor(exp)
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
### Create design matrix
###############################################################

design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
print(design)

###############################################################
### Estimate Dispersion
###############################################################

#y <- estimateDisp(y, design)

y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design, method="power")
y <- estimateGLMTagwiseDisp(y,design)

###############################################################
### Plot results
###############################################################


plotMDS(y)
plotBCV(y)


###############################################################
### Fit data
###############################################################

fit <- glmFit(y,design)

###############################################################
### Determine fold changes
###############################################################

mc <- makeContrasts(exp.r=WCFS1.glc-WCFS1.rib, levels=design)

fit <- glmLRT(fit, contrast=mc)

###############################################################
### Print top tags
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
###############################################################

DE_Genes<- fit$table[which(fit$table$PValue < 0.00005), names(fit$table) %in% c("logFC","LogCPM","LR","PValue")]

###############################################################
### joint de annotatie file met de DE_genes
###############################################################
library(dplyr)
names <- row.names(DE_Genes)
names <- data.frame(genes = names)
comp <- semi_join(WCFS1_anno,names, by = c("ORF" = "genes"))

sub <- select(filter(comp, comp$EC != ""),c("EC"))
comp2 <- semi_join(WCFS1_anno,sub, by = c("EC" = "EC"))

write.table(sub, file = "sub.txt")
write.csv(comp2, file = "MyData.csv",row.names=FALSE)
