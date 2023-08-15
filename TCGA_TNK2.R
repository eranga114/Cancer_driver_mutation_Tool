
#install.packages('BiocManager')
#library('BiocManager')
#BiocManager::install('TCGAbiolinks')
#BiocManager::install('maftools')
#BiocManager::install('pheatmap')
#BiocManager::install('SummarizedExperiment')


library(TCGAbiolinks)
library(tidyverse)
library(maftools)
library(pheatmap)
library(SummarizedExperiment)
library(biomaRt)


# Load TCGAbiolinks package
library(TCGAbiolinks)

# get a list of projects
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-COAD')


# build a query to retrieve gene expression data ------------
query_TCGA_test <- GDCquery(project = 'TCGA-COAD',
                            data.category = 'Transcriptome Profiling',
                            experimental.strategy = 'RNA-Seq',
                            workflow.type = 'STAR - Counts',
                            access = 'open')

getResults(query_TCGA_test)

# download data - GDCdownload
GDCdownload(query_TCGA_test)

# prepare data
tcga_brca_data <- GDCprepare(query_TCGA_test, summarizedExperiment = TRUE)
brca_matrix <- assay(tcga_brca_data, 'tpm_unstrand')

#convert into tibble
COAD_df = as_tibble(brca_matrix %>%as.data.frame()%>% rownames_to_column(var = "Gene"))

#Save as csv
write_csv(COAD_df, "COAD_df.csv", col_names = T)

# get a list of projects
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-STAD')


# build a query to retrieve gene expression data ------------
query_TCGA_test <- GDCquery(project = 'TCGA-STAD',
                            data.category = 'Transcriptome Profiling',
                            experimental.strategy = 'RNA-Seq',
                            workflow.type = 'STAR - Counts',
                            access = 'open')

getResults(query_TCGA_test)

# download data - GDCdownload
GDCdownload(query_TCGA_test)

# prepare data
tcga_brca_data <- GDCprepare(query_TCGA_test, summarizedExperiment = TRUE)
brca_matrix <- assay(tcga_brca_data, 'tpm_unstrand')

#convert into tibble
STAD_df = as_tibble(brca_matrix %>%as.data.frame()%>% rownames_to_column(var = "Gene"))

#Save as csv
write_csv(STAD_df, "TCGA-STAD_df", col_names = T)
