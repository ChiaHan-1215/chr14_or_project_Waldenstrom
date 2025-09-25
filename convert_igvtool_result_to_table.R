
###################################
# For the list of rsxxx.txt file  #
###################################

library(dplyr)

setwd('/Volumes/data/igvtool_count_GT_rstudio_biowulf_env/gt_res/')
ff <- list.files('.',pattern = '.txt')
df <- data.frame()

# need to modify column number when use! 

for (x in ff){
  tar <- read.delim(x)
  tar$snpid <- gsub('_.*','',x)
  tar$bam_ID <- gsub('.hg38','',tar$bam_ID)
  tar <- tar[,c(1,10,2,3,4,5,6)]
  df <- rbind(df,tar)
}



for (i in 1:nrow(df)){
  df$tot[i] <- sum(df[i,c(4:7)])
}

# remove row of total reads = 0 

for (x in 1:nrow(df)) {
  if (!is.na(df$tot[x]) && df$tot[x] == 0) {
    df <- df[-x, ]
  }
}

# 
# for (x in 1:nrow(df)){
#   for (y in 4:7){
#     if(as.numeric(df[x,y])/df$tot[x] > 0.1)
#     {df[x,y] <- gsub('','',names(df)[y])} else {
#       df[x,y] <- ""
#     }
#   }
# }


alleles <- c("A", "C", "G", "T")

for (a in alleles) {
  df[[paste0("called_", a)]] <- ifelse(df[[a]] / df$tot > 0.15, a, "")
}





#df <- df[,-8]
#df$snpid <- paste0(df$snpid,'_chr5:',df$pos)

df <- df[order(df$pos),]

f = df %>% group_by(pos) %>% tidyr::unite(col = target,c(9:12),na.rm = T,sep = '')
#f <- f[,-3]

# make A to AA 

for (x in 1:nrow(f)){
  if(f$target[x] %in% c('A','T','C','G'))
  {f$target[x] <- paste0(f$target[x]," ",f$target[x])}
  if (f$target[x] %in% c('AG','CT','AC')){
    f$target[x] <- gsub("([[:alpha:]])([[:alpha:]])", "\\1 \\2", f$target[x])}
}

# tansfrom dataset 
# df.wide <- tidyr::pivot_wider(f1, names_from = snpid, values_from = target)


df.wide <- f %>%
  ungroup() %>%
  tidyr::pivot_wider(
    id_cols = bam_ID,
    names_from = snpid,
    values_from = c(target, tot),
    names_glue = "{.value}_{snpid}"
  )


# rs149482608, A/G
# rs117972357, G/A
# rs117410836, T/C
# rs561799741, G/A
# rs534710671, C/A


# rs149482608 (A/G)
t1 <- f %>%
  filter(snpid == "rs149482608") %>%
  arrange(bam_ID) %>%
  transmute(
    bam_ID,
    rs149482608_A   = A,
    rs149482608_G   = G,
    total_read_rs149482608 = tot,
    rs149482608_gt = target
  )

t1 <- t1[,-1]


# rs117972357 (G/A)
t2 <- f %>%
  filter(snpid == "rs117972357") %>%
  arrange(bam_ID) %>%
  transmute(
    bam_ID,
    rs117972357_G   = G,
    rs117972357_A   = A,
    total_read_rs117972357 = tot,
    rs117972357_gt = target
  )

t2 <- t2[,-1]

# rs117410836 (T/C)
t3 <- f %>%
  filter(snpid == "rs117410836") %>%
  arrange(bam_ID) %>%
  transmute(
    bam_ID,
    rs117410836_T   = T,
    rs117410836_C   = C,
    total_read_rs117410836 = tot,
    rs117410836_gt = target
  )

t3 <- t3[,-1]

# rs561799741 (G/A)
t4 <- f %>%
  filter(snpid == "rs561799741") %>%
  arrange(bam_ID) %>%
  transmute(
    bam_ID,
    rs561799741_G   = G,
    rs561799741_A   = A,
    total_read_rs561799741 = tot,
    rs561799741_gt = target
  )

t4 <- t4[,-1]

# rs534710671 (C/A)
t5 <- f %>%
  filter(snpid == "rs534710671") %>%
  arrange(bam_ID) %>%
  transmute(
    bam_ID,
    rs534710671_C   = C,
    rs534710671_A   = A,
    total_read_rs534710671 = tot,
    rs534710671_gt = target
  )

t5 <- t5[,-1]


final <- t1 %>%
  left_join(t2, by = "bam_ID") %>%
  left_join(t3, by = "bam_ID") %>%
  left_join(t4, by = "bam_ID") %>%
  left_join(t5, by = "bam_ID")


## SAVE 

write.table(df.wide,'/Volumes/ifs/DCEG/Branches/LTG/Prokunina/CCLE_Whole_Genome_Sequencing_hg38/project_Waldenstrom_chr14_BAM_slice/ccle_gt_result.csv',col.names = T,row.names = F,quote = F,sep = ',')

write.table(final,'/Volumes/ifs/DCEG/Branches/LTG/Prokunina/CCLE_Whole_Genome_Sequencing_hg38/project_Waldenstrom_chr14_BAM_slice/ccle_gt_result_with_readcount.csv',col.names = T,row.names = F,quote = F,sep = ',')

