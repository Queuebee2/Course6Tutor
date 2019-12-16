### HAN (c) Todt
# coauthors SannevaSta, Queuebee2

library(edgeR)

###############################################################
### read counts
###############################################################

<<<<<<< Updated upstream
cnts <- read.delim("WCFS1_cnts.txt", comment.char="#")
WCFS1_anno <- read.delim2("WCFS1_anno.txt")

=======
# used datafiles
data_filename <- "WCFS1_cnts.txt"
annotation_filename <- "WCFS1_anno.txt"

# script directory
script_directory <-  dirname(toString(rstudioapi::getActiveDocumentContext()[2]))

# data directory
sep <- .Platform$file.sep
data_directory <-   paste0(sep, "data", sep)
full_data_path <- paste0(script_directory,data_directory)

# read files
cnts <- read.delim(paste0(full_data_path, data_filename), comment.char="#")
WCFS1_anno <- read.delim2(paste0(full_data_path, annotation_filename))
>>>>>>> Stashed changes

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

DE_Genes<- fit$table[which(fit$table$PValue < 0.00005), 
           names(fit$table) %in% c("logFC","LogCPM","LR","PValue")]

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

