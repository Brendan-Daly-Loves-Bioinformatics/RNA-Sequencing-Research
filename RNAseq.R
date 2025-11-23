# Load Libraries
library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

# Set the working directory
setwd()

# Load the count data
count_data <- read.csv('count_matrix.csv', header = TRUE, row.names = 1)
colnames(count_data)
head(count_data)

# Load the sample information
sample_info <- read.csv('design.csv', header = TRUE, row.names = 1)
colnames(sample_info)
head(sample_info)

# Set factor levels
sample_info$Treatment <- factor(sample_info$Treatment)
sample_info$Sequencing <- factor(sample_info$Sequencing)

# Create a deseq object and import the count data
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = sample_info, design = -Sequencing + Treatment)

# Set the reference for the Treatment factor
dds$Treatment <- factor(dds$Treatment, levels = c("untreated", "treated"))

# Filter the genes
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

# Perform the statistical tests to identify the differentially expressed genes
dds <- DESeq(dds)
deseq_result <- results(dds)
deseq_result

# Change DESeq object into an R object (dataframe)
deseq_result <- as.data.frame(deseq_result)
class(deseq_result)

head(deseq_result)

# Order the result table by increasing p values
deseq_result_ordered <- deseq_result[order(deseq_result$pvalue),]
head(deseq_result_ordered)

# Make queries
deseq_result[]






