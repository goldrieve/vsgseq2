# vsgseq2

An updated pipeline for analyzing VSG-seq data. The original VSGSeq pipeline is described in this [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4514441/) and [repository](https://github.com/mugnierlab/VSGSeqPipeline).

## Installation and Environment Setup

We use Nextflow/Conda to install dependencies, please install Nextflow as described [here](https://nf-co.re/docs/usage/installation) then run the following command to ensure everything is working:

```
nextflow run goldrieve/vsgseq2 -r main --help
```

Next, run the pipeline on synthetic demo data to ensure all dependencies have been installed. These can be downloaded, along with an example samplesheet.

```
wget https://github.com/goldrieve/vsgseq2/raw/refs/heads/main/data/reads/reads.tar
tar -xvf reads.tar
```

Now edit the samples.csv to point to the explicit location of the demo reads you just downloaded and run nextflow

```
nextflow run goldrieve/vsgseq2 -r main --samplesheet samples.csv
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

## vsgseq2 structure
The dag below summarises each step of vsgseq2.

![dag](figures/vsgseq2.dag)

## Customising analysis
It is possible to run sections of vsgseq2 using the --mode flag. The default is to run the whole pipeline, but say you have assembled the transcripts during a first run and wish to change a single flag in the analysis section, you can feed in the pre-assembled transcripts and start the pipeline from the analysis section. 

The default tutorial data run, including all vsgseq2 steps:
```
nextflow run goldrieve/vsgseq2 -r main --outdir tutorial_results
```

This is the same as running:
```
nextflow run goldrieve/vsgseq2 -r main --outdir tutorial_results --mode full
```

If you want to run the assembly step alone:
```
nextflow run goldrieve/vsgseq2 -r main --outdir tutorial_results --mode assemble
```

To run the entire analysis steps, post assembly:
```
nextflow run goldrieve/vsgseq2 -r main --outdir tutorial_results --mode analyse
```

The analyse steps can be broken down even further.
To predict VSGs from assemblies:
```
nextflow run goldrieve/vsgseq2 -r main --outdir tutorial_results --mode predictvsgs
```

To quantify the expression of predicted VSGs:
```
nextflow run goldrieve/vsgseq2 -r main --outdir tutorial_results --mode quantify
```

## Edit the pipeline execution using the following flags
```
VSGSEQ2.nf: A pipeline for analysing VSGSeq data

Required arguments:

  --assemblies Location of assemblies
                [default: data/tutorial_assemblies/*_trinity.Trinity.fasta]
  --vsg_db    Location of VSGdb
                [default: data/blastdb/concatAnTattb427.fa]
  --notvsg_db Location of NOTVSGdb
                [default: data/blastdb/NOTvsgs.fa]
  --vsgome    Location of VSGome
                [default: data/blastdb/concatAnTattb427.fa]
  --full_vsg_db    Location of a database to add into the VSGome (such as data/blastdb/concatAnTattb427_full.fa).
                    Default will only include the assembled VSGome.
  --mode    The mode to run the pipeline in. Options are full, assemble, predictvsgs, quantify, analyse.
                [default: full]


Optional arguments:

  --requestedcpus  Define number of cores VSGSeq2 will use.
                [default: 4]
  --cores  Define number of cores Trinity and other tools will use.
                [default: 4]
  --trinitymem  Define memory amount for Trinity in Gb.
                [default: 20 Gb]
  --cdslength    Define minimum CDS length (amino acids).
                [default: 300]
  --cdhitid       Define sequence identity threshold - how much the alignment has to match (0.0 - 1.0).
                [default: 0.98]
  --outdir        VSGSeq output directory.
                [default: results]
  --samplesheet  Define the path to the samplesheet.
                [default: samplesheet.csv]
  --help         Print this message.
  ```
