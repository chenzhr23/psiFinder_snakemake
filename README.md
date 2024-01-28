# psiFinder_snakemake

## Automated snakemake pipline with psiFinder quick start QT widget for Ψ-sites identification/annotation/target prediction

![quick_start](quick_start.png)

## Contents
- [Pre-installation](#pre-installation)
- [Input data and required files](#input-data-and-required-files)
- [Usage](#Usage)
  - [Set configuration](#set-configuration)
  - [Run psiFinder snakemake](#run-psiFinder-snakemake)

### Pre-installation
**psiFinder_snakemake** requires **snakemake/seqtk/cutadapt/STAR/bedtools/gtfToGenePred and several perl/python/R packages** pre-installation and predominantly used in unix-based operating systems. Therefore, for the usability of **psiFinder_snakemake**, we recommend running all the tools and scripts in WSL2 (WSL2 installation guide: https://pureinfotech.com/install-windows-subsystem-linux-2-windows-10/) or unix-based system with R and python.

Required perl modules:
```perl
use Getopt::Long;
use Bio::Perl;
use Bio::SeqIO;
```

Required python modules:
```python
import re
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
from sklearn.feature_extraction.text import CountVectorizer
import argparse
```

Required R packages:
```R
#use pacman to install packages in batch
install.packages("pacman")
library(pacman)

#load and install required R packages 
p_load("optparse","devtools","caTools","neuralnet","NeuralNetTools","dplyr","stringr","gridExtra","cowplot","pROC","mccr","ggplot2","ggpol","ggpubr","RColorBrewer","openxlsx","reshape2","factoextra","bedr","scales","e1071","tidyr")
```

### Input data and required files
Test data: 
- Input (CMC-control) file for Ψ-sites identification:**A1.cutadapt.extendedFrags.collapse.cutBarcodes.fa.gz** (download from: https://mega.nz/folder/oaUmhK7I#gSuYH4HW7OhL5qEbgmw0fw)
- Treat (CMC-treated) file for Ψ-sites identification:**A2.cutadapt.extendedFrags.collapse.cutBarcodes.fa.gz** (download from: https://mega.nz/folder/oaUmhK7I#gSuYH4HW7OhL5qEbgmw0fw)
- Genome files for STAR alignment:**hg38.fa**, **hg38.fa.fai** (download from https://mega.nz/folder/oaUmhK7I#gSuYH4HW7OhL5qEbgmw0fw and deposit it in  ./snakemake/genome)
- Annotation file for Ψ-sites identification: **hg38.gencode.v30.tRNA.refseqNcRNA.geneAnno.bed12** (download from https://mega.nz/folder/oaUmhK7I#gSuYH4HW7OhL5qEbgmw0fw and deposit it in  ./snakemake/script)
- Annotation file for Ψ-sites annotation: **hg38.genecode.v30.tRNA.snoRNA.miRNA.rmsk.exonFeatures.bed6** (download from https://mega.nz/folder/oaUmhK7I#gSuYH4HW7OhL5qEbgmw0fw and deposit it in  ./snakemake/script/annotation)
- Annotation file for Ψ-sites metagene: **gencode.v32.chr_patch_hapl_scaff.annotation.gtf** (download from https://mega.nz/folder/oaUmhK7I#gSuYH4HW7OhL5qEbgmw0fw and deposit it in  ./snakemake/script/metagene)

### Usage

#### Set configuration

![quickstart_config](quickstart_config.png)

#### Run psiFinder snakemake
```shell
bash run_psiFinder_snakemake.sh
```
