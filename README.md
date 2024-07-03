# vsgseq2
An updated  pipeline for analysing VSG-seq data. The original VSGSeq pipeline is described in this [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4514441/) and [repository](https://github.com/mugnierlab/VSGSeqPipeline).

## Installation and environment setup 

conda is used to set up the environment, please ensure you are using a recent conda version (anaconda 2024.02 and conda 24.3.0 work well).

```
git clone https://github.com/goldrieve/vsgseq2.git
cd vsgseq2
conda env create --file vsgseq2.yml -n vsgseq2
conda activate vsgseq2
```

## Quick start 
vsgseq2 is implemented using Nextflow, which is installed as part of the vsgseq2.yml.
To test the installation, use synthetic Illumina data for six samples which are stored in data/reads:
- paired reads:2 early, 2 late.
- single end reads: 1 early and 1 late.

To run vsgseq2 on the tutorial data, simply enter

```
nextflow run main.nf --outdir tutorial_results
```

This will create the directory __tutorial_results__ which will contain 4 subdirectories

1) VSGs 
- VSGs predicted for each sample (e.g. 1_VSGs.fasta).  
- concatenated list of all assembled VSGs (concatenated_vsgs.fasta). 
- final VSG database, after removing duplicate VSGs with cd-hit (VSGome.fasta).

2) assemblies 
- Trinity assembly for each sample.

3) summary 
- salmon alignment information (multiqc_report.html).
- quantification summary for each sample (tpm.csv).
- predicted VSG count for each sample (vsg_count.csv).

4) trimmed_reads 
- trimmed reads for each sample.

To visualise the expression data and number of assembled VSGs, use the R script bin/plot_script.R
Running the code will produce the figure below
![tutorial_figure](figures/tutorial_summary.png)

You can alter the run parameters using the following flags. A detailed DAG highlights where these flags will take effect - update!

## Input Files

main.nf takes both single and paired-end sequencing reads, in FASTQ format. Define the location of the reads in the samplesheet.csv and start the run!