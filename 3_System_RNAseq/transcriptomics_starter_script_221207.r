##Starterscript  

##################################################
##### Installation of libraries ##################
##################################################

##Install BiocManager
##For the prompts use "a" or "yes" (whatever fits)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

#Part 2
##Install Packages
##For the prompts use "a" or "yes" (whatever fits)
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
BiocManager::install("AnnotationDbi")
#BiocManager::install(organism, character.only = TRUE)
install.packages("ggridges")

#load libraries
# library(organism, character.only = TRUE)
library(clusterProfiler)
library("org.Hs.eg.db")
library(enrichplot)
library(GenomeInfoDb)
library(org.Hs.eg.db)
library(TBSignatureProfiler)
library(SummarizedExperiment)
library(pathview)
library(ggridges)
library(AnnotationDbi)

#Part 1
##Install Packages
##For the prompts use "a" or "yes" (whatever fits)
BiocManager::install("TBSignatureProfiler")
install.packages("tidyverse")
install.packages("matrixStats")
install.packages("limma")
#Load libraries
library(matrixStats)
library(tidyverse)
library(limma)

###################################################
##### Part 1: RNAseq data manual ##################
###################################################

##Load raw data data
hivtb_data <- TB_hiv

##Extract gene names and counts
counts = hivtb_data@assays@data@listData[["counts"]] # RO1 samples are TB positives

##################################################
##### Your code or from the script ###############
##################################################

boxplot(counts) # explora
logcounts <- log10(counts)
boxplot(logcounts)

mean_logcounts <- rowMeans(logcounts)
sd_logcounts <- rowSds(logcounts)

frame <- as.data.frame(cbind(mean_logcounts, sd_logcounts))

ggplot(data = frame, mapping = aes(x = mean_logcounts, y = sd_logcounts)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_bw()

plot(mean_logcounts, sd_logcounts)

q25 = quantile(mean_logcounts, 0.25)
which2save <- which(mean_logcounts > q25)
q25logcounts <- logcounts[which2save, ]

newdata <- q25logcounts[apply(q25logcounts, 1, IQR) > 1, ]
nrow(newdata)

# Data Exploration

tdata <- t(newdata)

pca <- prcomp(tdata, scale=T)
summary(pca)

plot(pca$x, type="n")
text(pca$x, rownames(pca), cex = 0.5)

# code the datapoints according to condtion TB or not
conditions = c(rep("TB", 16), rep("NOTB", 15))
plot(pca$x, type="n")
text(pca$x, labels = conditions, cex = 0.5)

# dendogram of correlation
pearsonCorr <- as.dist(1 - cor(newdata))
hC <- hclust(pearsonCorr)
plot(hC, labels = conditions)

heatmap(newdata)

# find the probes that differ significantly between TB status

condfactor <- factor(conditions) # factor levels
design <- model.matrix(~0+condfactor)
colnames(design) <- c("noTB", "TB")

## variance estimates with fit of data to design
fit <- limma::lmFit(newdata, design)
contrastmatrix <- limma::makeContrasts(TB - noTB, levels = design)

## p-values defined by contrast matrix
fit <- limma::contrasts.fit(fit, contrastmatrix)
ebayes <- limma::eBayes(fit)

hist(ebayes$p.value) # distribution of p-values

results <- limma::decideTests(ebayes)

## add gene symbols as row names
geneSymbol <- as.array(rownames((newdata)))
gs <- geneSymbol[c(which(results != 0))]
rownames(resData) <- gs

## filter for only results that passed the Test
resData <- newdata[results != 0, ]
pvalues <- ebayes$pvalue[results != 0]
resData <- cbind(resData, pvalues)

## add pvalues corrected
adj.pvalues <- p.adjust(ebayes$p.value, method = "BH")
adj.pvalues <- adj.pvalues[results != 0]
resData <- cbind(resData, adj.pvalues)

## output
write.table(resData, "most_regulated.txt", sep = "\t")



################################################################################
##### Part 2: Analysis with predefined pipelines from R packages ###############
################################################################################

##Alternative analysis with edgeR

##Load data again
hivtb_data = TB_hiv
##Extract counts per million 
hivtb_data = mkAssay(hivtb_data, log = TRUE, counts_to_CPM = TRUE)

##Look for genes with at least 3 samples with counts above 0.5 count per million
cpmdat <- cbind(hivtb_data@NAMES,hivtb_data@assays@data@listData[["counts_cpm"]]) %>% as.data.frame() %>% rowwise() %>% 
  mutate(tbsum = sum(across(R01_1:R01_9, ~sum(.x>0.5)))) %>% mutate(hivsum = sum(across(R02_17:R02_31, ~sum(.x>0.5)))) %>% 
  filter(hivsum > 2 & tbsum > 2) 

##extract raw counts for each genes and filter those out with less than 3 genes with 0.5 counts per million, 
##transform all samples columns into numeric variables
countdat <- cbind(hivtb_data@NAMES,hivtb_data@assays@data@listData[["counts"]]) %>% as.data.frame() %>% filter(V1 %in% cpmdat$V1) %>% 
  mutate(across(R01_1:R02_31, ~as.numeric(.x)))

##create edgeR object (group "1" = TB, group "0" = noTB)
egr = edgeR::DGEList(counts = as.matrix(countdat[,2:32]), group = c(rep(1,16),rep(0,15)))

##normalize fold changes in gene expression between samples
egr = edgeR::calcNormFactors(egr)

##Dispersion / distribution of gene counts
egr = edgeR::estimateDisp(egr)

##Common dispersion is assuming equal mean and standard-deviation among genes
##Tagwise dispersion assumes gene-specific dispersions (or of genes with a similar distribution)
edgeR::plotBCV(egr)

##Test gene expression differences between TB and noTB 
##The test is based on a negative binomial distribution 
et = edgeR::exactTest(egr)

et$table %>% 
  mutate(logp = -log10(PValue)) %>%
  ggplot(., aes(logFC,logp )) +
  geom_point(aes(color = logp), size = 2/5) +
  xlab(expression("log2 fold change")) + 
  ylab(expression("-log10 pvalue")) +
  scale_color_viridis_c()


####################################################
########### Gene enrichment analysis ###############
####################################################

#We need to specify our organism of interest, so humans.
organism = "org.Hs.eg.db"

##we need log2foldchange from the previous analysis and gene names
genes = cbind(et$table$logFC,cpmdat$V1) %>% as.data.frame()

#The gene names, as they are right now are in the "Gene card symbol" format.
#For the analysis we should changed them to the ensembl coding.
hs = org.Hs.eg.db
genenames = AnnotationDbi::select(hs, 
                   keys = genes$V2,
                   columns = c("ENSEMBL", "SYMBOL"),
                   keytype = "SYMBOL",
                   multiVals = "First")

#Some genes could not be recognized and are NA. Other have multiple ensembl ids.
#For simplicity remove the NAs and duplicates (take the first).
genes = genes %>% merge(.,genenames, by.x = "V2", by.y = "SYMBOL") %>% filter(ENSEMBL != "<NA>") %>% filter(!duplicated(ENSEMBL))

#Now extract the log2fold changes into a single vector, 
#name each value with the corresponding ensembl gene name, and sort the values decreasingly.
genenrich = as.numeric(genes$V1)
names(genenrich) = genes$ENSEMBL

#Sort the list in decreasing order 
genenrich = sort(genenrich, decreasing = TRUE)

#Now we can run the gene enrichment analysis
gse = gseGO(geneList=genenrich, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "fdr",
             eps = 0)

##Visualization of enrichment
require(DOSE)
##dot plot
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
##ridgeplot
ridgeplot(gse) + labs(x = "enrichment distribution")

####################################################
########### Bonus: KEGG Pathway analysis ###########
####################################################

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids=bitr(names(genenrich), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
#Remove duplicate ids again
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

#Create a new vector which has only the genes which were successfully mapped using the bitr function above
kegg_gene_list = genenrich[names(genenrich) %in% dedup_ids$ENSEMBL] 

#Name vector with ENTREZ ids
names(kegg_gene_list) = dedup_ids$ENTREZID

#Omit any NA values 
kegg_gene_list=na.omit(kegg_gene_list)

#Sort the list in decreasing order again
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

##Humans organism again
kegg_organism = "hsa"

##kegg pathway analysis
kk2 = gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid",
               eps = 0)
id_list <- list(as.data.frame(kk2)$core_enrichment)


#Produce the KEGG plott (look for the png file written by the below function on your computer)
#Choose "pathway.id" to the id of your interest. Look at kk2 to extract relevant ones.
dme = pathview(gene.data=kegg_gene_list, pathway.id= "05152", species = kegg_organism)
