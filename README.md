# chr14_or_project_Waldenstrom

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
