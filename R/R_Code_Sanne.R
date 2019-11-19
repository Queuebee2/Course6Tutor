# maakt twee groepen dus glucose en ribose
group <- c(1,1,2,2)
# maakt een bewerkbare lijst van de dataset
y <- DGEList(counts=x, group=group)

# filtert genen die bijna niet tot expressie komen eruit
keep <- filterByExpr(y)
# hou alleen degene die dus bruikbaar zijn
y <- y[keep, , keep.lib.sizes=FALSE]

# normalizatie
y <- calcNormFactors(y)

# maakt een plot (The function plotMDS produces a plot in which distances between samples
#correspond to leading biological coefficient of variation (BCV) between those samples)
plotMDS(y)

# maakt de data voor dispersie plot  maar er moet nog een design bij  maar snap niet hoe ik dat moet maken dus deze klopt wss niet
y <- estimateDisp(y, robust = TRUE)

# plot of dispersion
plotBCV(y)