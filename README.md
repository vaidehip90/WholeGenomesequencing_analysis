# Whole Genome sequencing analysis
Whole genome sequencing is method for analyzing the entire genome. Whole genome sequencing is used for variety of purpose such as identifying genetic disorder,characterizing the mutations and also foodborne pathogen identification by FDA. In this project, I used three datset from 1000 genome project phase 3 to identify the somatic variants present in the dataset. The aim of this study is to identify specific somatic variation present with all information such as genes, chrosome location, variants,SNPs or del var

The Whole genome sequencing pipeline includes quality check by fastqc, alignment to reference genome by BWA-MEM,Mark Duplicate tool to remove the duplicate reads,gatk Mutect2 variant caller to call somatic variants,gatk Funcotator to annotate the variants.
.In this workflow,the bioinformatics tools were installed using conda.

# Whole genome sequencing workflow

> STEP1: Downloading the data and reference genome

> STEP2: Fastqc

> STEP3: BWA MEM

> STEP4 :Mark duplicate

> STEP5: Mutect2

> STEP6: Gatk Funcotator 
