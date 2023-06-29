rm(list = ls())

path_to_files = "6_16SRNA/fastqs/" 
list.files(path_to_files) # lists all files

# install.packages(tidyverse)
BiocManager::install("ShortRead")
library(tidyverse)
library(ShortRead)
fastq = readFastq(paste0(path_to_files,list.files(path_to_files)[1]))

fastq@sread # displays reads

# determine bacteria origin

BiocManager::install("dada2")
BiocManager::install("phyloseq")
library(dada2)
library(phyloseq)

# forward reads
fnFs = sort(list.files(path_to_files, pattern="R1_001.fastq", full.names = TRUE))

# reverse reads
fnRs = sort(list.files(path_to_files, pattern="R2_001.fastq", full.names = TRUE))

# extracting the sample names
sample_names = gsub("JG-JG-","JG-",fnFs)
sample_names = sapply(strsplit(sample_names, "-"), `[`, 2)
sample_names = sapply(strsplit(sample_names, "_"), `[`, 1)

meta_dat = read.csv("6_16Srna/meta_tb.csv") #%>% dplyr::filter(id != "")

# specifiy filtered datapaths for later

path_to_files_fil = "6_16Srna"
filtRs = file.path(path_to_files_fil, "filtered", paste0(sample_names, " R_filt.fastq.gz"))
names(filtRs) = sample_names
filtFs = file.path(path_to_files_fil, "filtered", paste0(sample_names, " F_filt.fastq.gz"))
names(filtFs) = sample_names

# select a subset of all samples

metadata_use = meta_dat %>% filter(Lab.ID %in% c("6916007",
                                                 "6916008",
                                                 "6916009",
                                                 "6916010",
                                                 "6917492", 
                                                 "6917490",
                                                 "1515745",
                                                 "6917657"))

fnFs = fnFs[names(filtFs) %in% metadata_use$Lab.ID]
filtFs = filtFs[names(filtFs) %in% metadata_use$Lab.ID]
fnRs =fnRs[names(filtRs) %in% metadata_use$Lab.ID]
filtRs = filtRs[names(filtRs) %in% metadata_use$Lab.ID]

# Quality/Q-scores of reads, read position vs Q

# plotQualityProfile(fnFs)
# plotQualityProfile(fnRs)
# beginning has very high Q score compared to end
# 0, 0 and 30

out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen= c(200,200), maxN=0, maxEE=2, truncQ=20, 
                    rm.phix=TRUE, compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

errR = learnErrors(filtRs, multithread=FALSE)
errF = learnErrors(filtFs, multithread=FALSE)

plotErrors(errR, nominalQ=T)
plotErrors(errR, nominalQ=T)

# extract sample composition

dadaRs = dada(filtRs, err = errR, multithread=F)
dadaFs = dada(filtFs, err = errF, multithread=F)

# final cleaning and data formating

mergers = mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# rows as samples, columns as unique sequence 

seqtab = makeSequenceTable(mergers)

#  recognizes chimeric sequences and removes those

seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=F, verbose=TRUE)

df <- errR$trans[, 22:38]
x <- colSums(df)

for (i in 1:ncol(df)){
  df[, i] <- df[, i] / x[i]
}

df_base_error <- rowSums(df) 

matrix <- matrix(data = df_base_error, nrow = 4, ncol = 4,byrow = T,
                 dimnames = list(c("A", "C", "G", "T"), c("A", "C", "G", "T")))

for (i in 1:4){
  matrix[i, i] <- 0
}

base::heatmap(matrix, Rowv = NA, Colv = NA)

legend(x = "right", legend = c("max", "min", "correct"), fill = heat.colors(4))
