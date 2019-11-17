# todo

- [ ] install anaconda jupyter notebook environment with R
- [ ] do stuff
- [ ] do ore stuff


1. Installaties
- [x] Install R
- [x] install bioconductor
-  [x] install specific packages from bioconductor
  - [x] bioBase
  - [x] edgeR
  - [x] DEseq
  - [x] DEseq2
2. User manuals
Zoals ik het nu begrijp, zijn user manuals verstopt. Daarvoor moeten we met de openVignette() functie, van bioBase, een package's "Vignette"(?) openen. Bijvoorbeeld `bioBase::openVignette("edgeR")`. Daarvoor moet de package wel eerst geladen zijn. Dat doe je met `library(<packagename>)` Wat er precies opent begrijp ik niet echt. Voor EdgeR komt er bijvoorbeeld een korte pagina met daarin nog een command om uit te voeren om bij de 'echte' manual te komen.
3. data bekijken
laad de data in met `dat <- read.delim(file="RNA-Seq-counts.txt", row.names=1)
`
  - [dit is interessant](https://www.biostars.org/p/130044/)
  -

# bioBase en andere packages installeren

RStudio moet geopend zijn met administratorrechten. Daarna kun je met de volgende commands bioBase installeren (source::bioconductor/biobase)
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("bioBase")
```
Andere packages zijn makkelijk te installeren door `"bioBase"` te vervangen met de packagenaam (let op case sensitivity) zoals `"edgeR"`,`"DESeq"`,`"DESeq2"`... etc.

Na het installeren van een package (en minimaal bioBase), k
```
library(bioBase)              # this is needed to use the openVignette
library(edgeR)
library(DESeq)
library(DEseq2)

openVignette("<packagename>") # or instead of loading bioBase entirely
                              # bioBase::openVignette("<packagename>")
```


### normalisatiefuncties per package  
##### EdgeR
Normalisation in edgeR is concerned with relative changes[^1]




[^1]: edgeRManual
