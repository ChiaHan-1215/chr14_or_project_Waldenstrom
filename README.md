# chr14_for_project_Waldenstrom

### Goal: extract Genotype files from CCLE, GTEx, CARD for finding rare variants for Waldenstrom project

the rare variants are located in chr14

```
RS_Number      Position (hg38_HC)     Allele Frequency
rs149482608    chr14:95525224 A=0.998, G=0.002
rs117972357    chr14:95577209 G=0.996, A=0.004
rs117410836    chr14:95585637 T=0.971, C=0.029
rs561799741    chr14:95763084 G=0.998, A=0.002
rs534710671    chr14:95845144 A=0.998, C=0.002
```

### Code: 

- Extract BAM slice from CCLE using samtools in Terra.bio 

**********

### Future direction:

- Extracted some cell of interest around target variants 
https://genome.ucsc.edu/s/gahanleeo/hg19_chr14_project_Waldenstrom

- Search B cell lymphoma HiCHIP data, list link of interest:

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136090

#### sc-RNA/atac for Waldenstrom patient
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296167

- This dataset have some of interest dataset, trying to get bulked snATAC-seq and output as bigwig for display
- `ArchR` package has function to do that, code in biowulf.

```R

# Goal: Using ArchR to convert sample sc-ATAC to Bigwig and display in UCSC genome broswer
# Date: 10062025

library(dplyr)
library(ArchR)
# ArchR::installExtraPackages()

setwd('GSE296167_snATAC_track_convert_to_bigwig/')

addArchRThreads(threads = 4) 
addArchRGenome("hg38")

inputFiles <- list.files('.')
names(inputFiles) <- c("S1","S2","S3","S4")


# carte ArchR file 
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

#  Create an ArchRProject (choose an output folder)

proj <- ArchRProject(ArrowFiles, outputDirectory = "ArchROut", copyArrows = TRUE)

# Export one BigWig per Sample (pseudo-bulk coverage)
bw <- getGroupBW(
  ArchRProj = proj,
  groupBy   = "Sample",
  normMethod= "ReadsInTSS",     # recommended
  tileSize  = 100,              # bin size; smaller = smoother but larger files
  maxCells  = 1000              # subsample cells per group; increase if desired
)

# check, the outout will located in the work dictionary as folder
bw

```



