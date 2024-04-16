# vsgseq2
An updated  pipeline, implemented using nextflow, for analysing VSG-seq data. Originally descrived in this [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4514441/).

## Basic Usage 
```
git pull https://github.com/goldrieve/vsgseq2
cd vsgseq2
nextflow run vsgseq2.nf
```

You can adjust various parameters for your run. Read about all of the options using `nextflow run vsgseq2.nf --help`.


|====================|
| V S G S E Q 2 - N F|
|====================|

VSGSEQ2.nf: A pipeline for analysing VSGSeq data

Required arguments:

  --reads Location of reads, if not in reads dir
                [default: /Users/goldriev/pkgs/vsgseq2/data/reads/*{1,2}.fq.gz]
  --vsg_db    Location of VSGdb
                [default: concatAnTattb427.fa]
  --NOTVSG_db Location of NOTVSGdb
                [default: NOTvsgs.fa]
  --bait Location of reference genome to be used as bait for salmon quantification
                [default: TriTrypDB-67_TbruceiEATRO1125_Genome]

Optional arguments:

  --requestedcpus  Define number of cores VSGSeq2 will use.
                [default: 4]
  --cores  Define number of cores Trinity will use.
                [default: 4]
  --trinitymem    Define mem Trinity will use.
                [default: 20]
  --cdslength    Define minimium CDS length (amino acids).
                [default: 100]
  --cdhitid       Define sequence identiy threshold - how much the alignment has to match (0.0 - 1.0).
                [default: 0.98]
  --outdir        VSGSeq outdir.
                [default: results]

## Dependencies

Dependencies should be installed using anaconda


## Input Files

vsgseq2 has been updated to allow to submissiong of paired-end sequencing reads, in FASTQ format. Place the FASTQ files in the directory data/reads and the pipeline will do the rest. If you want to run the pipeline on a subset of these reads, add the flag --reads 'dir/*{1,2}.fq.gz' specifiying the location and names of files you would like to analyse.

## BLAST Databases for Identifying VSGs

You can use any reference you want to identify VSGs. We have a few options available in data/blastdb/.

We provide the VSG database:
* Combined database of BOTH Lister427 and EATRO1125 (`concatAntattb427`)

To create your own, place the fasta file in the data/blastdb/ directory, index it with blast and provide the following flag upon job submission --vsg_db new_db

There is one 'NonVSG' database (`NOTvsgs`). This database has been cobbled together after multiple iterations of assembling expressed VSGs, inspecting them by hand, and identifying common false positives (e.g., certain ESAGs assemble frequently and this will filter those out). If you run the pipeline using this filter (the default), you'll need this database available for BLAST.  

## Output Files

Output files are saved in one folder '--outdir results'. A summary file shows the expression of each VSG in each sample, both in terms of TPM and percentage of the population (TPM for that VSG/total TPM).

