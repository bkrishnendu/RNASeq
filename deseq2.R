#install.packages("usethis")
#usethis::edit_r_environ()
#usethis::edit_r_profile()
#Add to your .Renviron the following string:
#CURL_SSL_BACKEND=openssl
#Add to your .Rprofile the following string:
#options(download.file.method="libcurl", url.method="libcurl")

#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
BiocManager::version()
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
BiocManager::install("ComplexHeatmap")
BiocManager::install("GEOquery",force = TRUE)

library("DESeq2")
library("edgeR")
library("dplyr")
library("tidyverse")
library("RColorBrewer")
library("ComplexHeatmap")



##Defining the data directory
data_dir<-"F:/RNA/GSE225620/"
out_dir <-"F:/RNA/GSE225620/Output/" 


##Pre-processing the data frame for edgeR
raw_mat<-read.csv("GSE225620/GSE225620_GeneLevel_Raw_data.csv", header=TRUE, row.names = 1)

raw_mat<-raw_mat %>% remove_rownames() # remove the ensembl_id column

if(length(which(duplicated(raw_mat$gene_symbol)))>=1){
  print("Duplicated rows found and removing them")
  raw_mat<-distinct(raw_mat,gene_symbol, .keep_all = TRUE)
  print("Duplicated rows removed")
} else{
  print("No duplicated rows found")
}
raw_mat<-raw_mat %>% remove_rownames %>% column_to_rownames(var="gene_symbol")


## Read the clinical_data
clinical_data<- read.csv("GSE225620/SraRunTable.csv", header=TRUE, row.names =1)

#Remove the normal samples
clinical_data <- clinical_data %>% filter(!timepoint %in% "Normal")

table(clinical_data$timepoint)



## Take the common sample names match in both clinical and raw expression data
commonIDs <- intersect(rownames(clinical_data), colnames(raw_mat))
raw_mat = raw_mat[,commonIDs]
clinical_data = clinical_data[commonIDs, ]

## Check that sample names match in both files
all(colnames(raw_mat) %in% rownames(clinical_data))
all(colnames(raw_mat) == rownames(clinical_data))

## Create DESeq object
dds <- DESeqDataSetFromMatrix(raw_mat, colData=clinical_data, ~timepoint)

##For comparison below, we start with minimal filtering with edgeR.
#This basically removes genes with very low counts across most samples.
y <- DGEList(counts=counts(dds), genes = row.names(raw_mat), 
             group = factor(clinical_data$timepoint))
keep <- filterByExpr(y)
table(keep)
y <- y[keep,]
dds <- dds[keep,]

## Set the 'specimen_type' variable as a factor with defined levels
# This ensures that 'post-treament' is the first level and 'pre-treament' is the second level in the analysis
# Specifying levels helps control the order in which the groups are compared in differential expression analysis
dds$timepoint <- factor(dds$timepoint, levels = c("post-treament", "pre-treament"))

## Run the DESeq2 analysis pipeline on the DESeqDataSet object 'dds'
# This function performs the differential expression analysis
# It includes steps such as:
# 1. Estimating size factors for normalization
# 2. Estimating dispersion for each gene
# 3. Fitting a generalized linear model (GLM) for each gene based on the specified design
# 4. Performing statistical tests for differential expression
# The result is stored back in the 'dds' object, which now contains additional information
# such as normalized counts, dispersion estimates, and results of the statistical tests
dds <- DESeq(dds)

resultsNames(dds)

## Plot dispersions : Plotting the dispersion estimates is a useful diagnostic. 

#The dispersion plot below is typical, with the final estimates shrunk from the gene-wise estimates towards the fitted estimates.

#The dispersions plot is a key diagnostic plot that visualizes the relationship between gene expression mean and dispersion (or variability) across samples. 
#This plot is important because DESeq2 models gene expression counts using a negative binomial distribution, where the dispersion parameter describes how much the counts deviate from the average (mean) for each gene.

#1.  Dispersion: Dispersion measures the variability of gene expression across replicates for a given gene. Higher dispersion means that the expression of a gene varies more between samples. Genes with low mean expression often have higher dispersion, while highly expressed genes tend to have lower dispersion.
#2.  Mean-Dispersion Relationship: The plot shows the estimated dispersions against the mean normalized counts for each gene:
  #- X-axis: Mean expression level (average normalized counts) for each gene.
  #- Y-axis: Estimated dispersion for each gene.
#3.  Fitted Line: DESeq2 fits a curve to model the relationship between mean expression and dispersion. The plot includes a red line representing this fitted model. The model accounts for the expected variability of genes based on their mean expression, allowing for more accurate differential expression analysis.
#4.  Gene-Specific Estimates: The black points in the plot represent the gene-specific dispersion estimates, which reflect the raw variability in the data before shrinkage (correction).
#5.  Shrunken Dispersion Estimates: The blue points represent the final, shrunken dispersion estimates after DESeq2 applies shrinkage to stabilize the dispersion estimates, particularly for lowly expressed genes. Shrinkage helps to prevent overestimation of variability for genes with low counts, which might otherwise lead to false positives in differential expression analysis.

plotDispEsts(dds, main="Dispersion plot")

## Perform variance stabilizing transformation
res <- results(dds,contrast=c("timepoint","post-treament","pre-treament"))
summary(res)

res <- as_tibble(rownames_to_column(as.data.frame(res), "gene"))

write.csv(res, file="GSE225620/deseq2_degs.csv", row.names = F)


## Save the normalised data
# Extract the normalized counts from the DESeqDataSet object 'dds'
# The 'counts' function retrieves the normalized counts when the argument 'normalized=T' is set
norm_mat <- counts(dds, normalized=TRUE)

write.csv(as.data.frame(norm_mat),file=paste0(out_dir,"deseq2_normalized_counts.csv"))

## Perform variance stabilizing transformation
# Perform variance stabilizing transformation (VST) on the DESeqDataSet object 'dds'
# This transformation stabilizes the variance across samples, making the data more suitable for clustering and visualization
# 'blind = FALSE' means that the transformation takes the experimental design into account
vsd <- vst(dds, blind = FALSE)

##Plot PCA
plotPCA(vsd, ntop=500, intgroup ="timepoint")



## Plot Heatmap of Sample Distances
# Calculate the distance matrix based on the transformed counts
# 'dist(t(assay(vsd)))' computes the pairwise distances between samples
# Here, we apply the dist function to the transpose of the transformed count matrix to get sample-to-sample distances

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)


# Create a data frame of annotations
ann <- data.frame(
  timepoint = colData(vsd)$timepoint,
  sex = colData(vsd)$sex,
  row.names = colnames(assay(vsd))
)


# Define colors for the annotations
timepoint_colors <- setNames(c("orange2", "darkblue"), unique(ann$timepoint))

sex_colors <- setNames(c("orange2", "darkblue"), unique(ann$sex))

# Create a heatmap annotation for the columns

colAnn <- HeatmapAnnotation(
  timepoint = ann$timepoint,
  sex = ann$sex, # Add the 'sex' annotation
  col = list(
    timepoint = timepoint_colors,
    sex = sex_colors
  )
)

# Create a heatmap annotation for the columns
rowAnn <- rowAnnotation(
  timepoint = ann$timepoint,
  sex = ann$sex,
  col = list(
    timepoint = timepoint_colors,
    sex = sex_colors
  )
)

# Define the color palette for the heatmap
myCol <- colorRampPalette(c("blue", "black","yellow"))(100)


# Plot the heatmap
Heatmap(
  sampleDistMatrix,
  name = "Distance",
  top_annotation = colAnn,  # Add the column annotation
  left_annotation = rowAnn,
  col = myCol ,             # Specify the color palette
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,       # Enable row clustering
  cluster_columns = TRUE,   # Enable column clustering"
  row_names_side = NULL
)

