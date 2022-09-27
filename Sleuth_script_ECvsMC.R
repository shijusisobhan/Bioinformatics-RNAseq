# sluth analysis Raw data


# Input --> kallisto abundance table
#sample    transcripts    scaled     reads_per_base      tpm


# Output: beta value--> log2 fold changes between conditions
#ens_gene: associated Ensembl gene ID
#ext_gene: associated gene symbol



# Clear the workspace

rm(list=ls())

# load the required libraries

library(sleuth)
library(biomaRt)
library(devtools)
library('ggplot2')


# change the directory: Please change the following lines depends on your file location

setwd('D:/Sleuth_analysis')

base_dir <- 'D:/Sleuth_analysis/Kallisto_out'


# Sample name
sample_id <- dir(file.path(base_dir))

# Location of alignment output for each sample
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))

# Read the sample-condtion information table (We have to made it based on experiment)
s2c <- read.table("SC_ECvsMC.txt", header = TRUE, stringsAsFactors=FALSE)

# Add a column that locate the directoies of each sample and condtion. This column must be labeled path
# The user should check whether or not the order is correct. 

s2c <- dplyr::mutate(s2c, path = kal_dirs)




#********************get the gene names using biomaRt************************************


# collect gene names  from Ensembl

#mart <- useEnsembl(biomart = "ensembl", dataset = "dmelanogaster_gene_ensembl")


mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "dmelanogaster_gene_ensembl", host="uswest.ensembl.org",ensemblRedirect = FALSE)


# Add genename into the sleuth table
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id","external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,ens_gene = ensembl_gene_id, ext_gene = external_gene_name)



so<- sleuth_prep(s2c, ~ condition, target_mapping = t2g,aggregation_column = 'ext_gene', gene_mode = TRUE, extra_bootstrap_summary=T,
                 read_bootstrap_tpm = TRUE)


#*************** fitting the model**********************************************

# fitting our data to distributions (Data--> Number of reads of transcripts)

# (2) estimate parameters for the sleuth response error measurement model

# read_count ~ condition 

so <- sleuth_fit(so, ~condition, "full")

# wald test

# We can use a Wald test to measure the effect of a particular condition against the control 
#condition for a variable.

# Wald tests calculate p values to estimate the significance of any difference between conditions.
# Small p value--> significant

# 1.How do we know what models are available for testing?
models(so)

# run Wald tests on each of the coefficient (beata)

beta='conditionMC'

so <- sleuth_wt(so, which_beta = beta)


# Fold change (b) value

# Sleuth fits log normal distributions using natural logs
#To get the raw effect sizes, take the exponents of the beta values;
#you could then calculate log2 effect sizes:

# Fold change (b) value

# Sleuth fits log normal distributions using natural logs
#To get the raw effect sizes, take the exponents of the beta values;
#you could then calculate log2 effect sizes:

sig_level=0.1
test_table <- sleuth_results(so, beta)
test_table$raw_b <- exp(test_table$b)
test_table$log2_b <- log2(test_table$raw_b)
test_table$neg_log10_qval<- -log10(test_table$qval)


# get the gene expression table

sleuth_gene_matrix <- sleuth_to_matrix(so, 'obs_norm', 'tpm')

DGE_table<-cbind(row.names(sleuth_gene_matrix), sleuth_gene_matrix)# Add row names as one colum (gene names

p_val_sort<-(test_table[order(test_table$target_id),]) # sort according to gene name


p_val_reduce<-p_val_sort[c("target_id","pval","qval", "log2_b")] # extract only genename,p and q and log2()


DGE_all<-merge(DGE_table, p_val_reduce, by.x=1 , by.y="target_id") # megrge with gene expression and p, q values

DGE_all<-na.omit(DGE_all)

#DGE_all<-DGE_all[which(DGE_all$qval<0.1 & abs(DGE_all$log2_b)>1 ), ]


DGE_all<-(DGE_all[order(-DGE_all$log2_b),])

write.csv(DGE_all, 'DGE_table_ECvsMC.csv')



test_table$diffexpressed <- "Not sig"
# if log2Foldchange > 0 and qvalue < 0.1, set as "UP" 
test_table$diffexpressed[test_table$log2_b > 0 & test_table$qval < sig_level] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
test_table$diffexpressed[test_table$log2_b < 0 & test_table$qval < sig_level] <- "DOWN"


#****************************************************************************************

N_significant<-length(test_table$diffexpressed[test_table$diffexpressed !="Not sig"])
N_UP<-length(test_table$diffexpressed[test_table$diffexpressed =="UP"])
N_DOWN<-length(test_table$diffexpressed[test_table$diffexpressed =="DOWN"])

N_significant
N_UP
N_DOWN


#********************************************************************************

library('ggplot2')
ggplot(test_table) + geom_point(aes(x = log2_b, y = neg_log10_qval, col=diffexpressed))+
  geom_vline(xintercept=0, col="black",  linetype="dashed")+
  scale_color_manual(values=c("blue", "black", "red"))+
  labs(x = "log2(FC)", y="-log10(q)", colour="DEG")
