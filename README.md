# Bioinformatic steps in RNA sequencing (Tutorial)

Quantifying gene expression and identifying transcripts that are differentially expressed between two set of samples is an important approach in modern biotechnology. RNA-seq analysis based on next-generation sequencing (NGS) data is the gold standard for the analysis of gene expression at the level of the whole transcriptome.  RNA-seq involves isolation of total RNA from tissues or cells of interest followed by the construction of DNA libraries and sequencing of these libraries using a next-generation sequencing instrument (Fig 1). 

![Fig1](https://raw.githubusercontent.com/shijusisobhan/Bioinformatic-RNAseq/main/Figures/Fig1.jpg)
*Fig 1. RNA sequencing steps (From Griffith et al. (2015))*

In this tutorial we explain the different bioinformatic steps involved in the analysis of RNA seq data to find the differentially expressed genes in two distinct set of samples. Fig 2 shows the current work flow of bioinformatic analysis of RNA seq data. 

![Fig2](https://raw.githubusercontent.com/shijusisobhan/Bioinformatic-RNAseq/main/Figures/Fig2.jpg)
*Fig 2. Current RNA seq analysis work flow*

## Data selection for the analysis
Drosophila is a fruit fly belonging to the family Drosophilidae. Drosophila melanogaster, one species of drosophila, has been widely used as a model organism in developmental biology. A study has been conducted by Allada Lab to understand the molecular mechanisms by which discrete circadian clock neurons program a homeostatic sleep center [Full Text](https://www.biorxiv.org/content/10.1101/2021.10.22.465404v1.full). In normal cycles of day and night, Drosophila exhibit morning and evening peaks of activity, which are controlled by two different group of neurons in the brain. In this study RNA sequence has been done and collected the data from morning and evening cells. There are 12 samples in this study and all are deposited in the Gene Expression Omnibus (GEO) with an accession number GSE186076 [GEO accession](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186076). In this tutorial we are looking on only 6 samples Morning control (3 samples MC1, MC2, MC3) and Evening control (3 samples EC1, EC2, EC3).  We select these data set, because these are some special case which will explained in the proceeding sections. Our aim is to performing RNA seq analysis on these data and find out differentially expressed genes in Morning and evening cells. 

### I. Download the RNA seq raw data (fastq, fasta…) from the public data base.
There are many methods are existing to download raw data from the public data base. Here I am describing one of the easiest methods to download the data. Users can use any method to download the data. 

Note: Raw data are usually big in size (normally GBs). If the user is not interested in downloading the data or your internet connection has only limited bandwidth, you can skip this downloading step. Any way I am providing the kallisto output [here](https://github.com/shijusisobhan/Bioinformatic-RNAseq) and user can download it directly from the repository and perform sleuth analysis (Differential gene expression analysis).

Downloading steps
1. Go to the European Nucleotide Archive (ENA) [website](https://www.ebi.ac.uk/ena)
2. In the search tab on the top right, enter the GEO accession number  GSE186076

![Fig3](https://raw.githubusercontent.com/shijusisobhan/Bioinformatic-RNAseq/main/Figures/Fig3.jpg)
*Fig 3*

3. ENA provide following search results

![Fig4](https://raw.githubusercontent.com/shijusisobhan/Bioinformatic-RNAseq/main/Figures/Fig4.jpg)
*Fig 4*

 4. Click on the SRP341965 tab, you will be directed to following page
 
 ![Fig5](https://raw.githubusercontent.com/shijusisobhan/Bioinformatic-RNAseq/main/Figures/Fig5.jpg)
 *Fig 5*
 
 Now user can download all the files using Download All option or selected files by checking the box near the files. 
 
Note: In this example there are two files per each entry, for example SRR****_1.fastq.gz & SRR****_2.fastq.gz. These are paired end sequencing (pair-1, pair-2). Kaliisto analysis is different for paired end sequencing and single end sequencing. I will explain it well in the kallisto analysis section. 

In this tutorial we are only interested on EC and MC samples. When you looking at the read file window above, samples names are missing, which is difficult to understand the users. So if you need the sample names click on the ‘show column section’ tab and check the ‘sample_title’ box. Then slide the window, you can see the samples name on the right.
![Fig6](https://raw.githubusercontent.com/shijusisobhan/Bioinformatic-RNAseq/main/Figures/Fig6.jpg)
*Fig 6*

Now you can easily identify the samples.











### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [Basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/shijusisobhan/Bioinformatic-RNAseq/settings/pages). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.
