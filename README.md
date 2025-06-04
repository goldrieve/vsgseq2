# vsgseq2

An updated pipeline for analyzing VSG-seq data. The original VSGSeq pipeline is described in this [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4514441/) and [repository](https://github.com/mugnierlab/VSGSeqPipeline).

## Installation and Environment Setup

We use Nextflow to run vsgseq2 analysis. You can install Nextflow via the manual vsgseq2 installation method.

1) Manual installation with conda used to install dependencies:

```
git clone https://github.com/goldrieve/vsgseq2
cd vsgseq2/
conda env create -f vsgseq2.yml
nextflow run main.nf --help
```

2) Installation via Nextflow

Nextflow can directly install vsgseq2 with the help of Docker, Singularity or Conda to install dependencies. Without cloning the vsgseq2 github, run:

```
conda create --name nf-env bioconda::nextflow
nextflow run goldrieve/vsgseq2 -r main --help
```

Either way, you can test the installation with synthetic VSGSeq data:

```
wget https://github.com/goldrieve/vsgseq2/raw/refs/heads/main/data/reads.tar.gz
tar -xzf reads.tar.gz
cd reads
```

Once you have navigated to the synthetic reads directory, run vsgseq2 on the data with the provided samplesheet. Here we use Docker, Conda and Singularity are also options which we describe below.

```
nextflow run goldrieve/vsgseq2 -r main -with-docker --samplesheet samples.csv
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

```mermaid
---
config:
  themeVariables:
    fontSize: 30px
---
flowchart TD
 subgraph VSGSeq["."]
        C2{"_Trim reads?_"}
        D2["Trim Galore & Cutadapt"]
        E2["Use raw FASTQ files"]
        F2["Trinity assembly"]
        G2["Find ORFs in contigs"]
        H2["Concatenate ORFs"]
        I2["Cluster with cd-hit-est"]
        J2["BLAST merged ORFs vs VSG DB"]
        M2["Map reads with Bowtie"]
        N2["Quantify with MULTo"]
        P2["Write final results"]
  end
 subgraph s1["."]
        C1{"_Trim reads?_"}
        D1["Trim Galore & Cutadapt"]
        E1["Use raw FASTQ files"]
        F1["Trinity assembly"]
        G1["Find ORFs in contigs"]
        H1["Concatenate ORFs"]
        I1["Cluster with cd-hit-est"]
        J1["BLAST merged ORFs vs VSG DB"]
        M1["Map reads with Bowtie"]
        N1["Quantify with MULTo"]
        P1["Write final results"]
        Q1["Files are analysed in parallel"]
  end
 subgraph vsgseq2["."]
        C3["Trim reads"]
        D3{"_Full or Analyse?_"}
        E3["Trinity assembly"]
        F3["Find ORFs in contigs"]
        G3["BLAST merged ORFs vs VSG DB"]
        H3["Concatenate VSGs"]
        I3["Cluster with cd-hit-est"]
        K3["Quantify with Salmon"]
        L3["Write final results"]
        M3["Multiple files are analysed in parallel"]
  end
    File1["File1"] --> Q2["Files are analysed independently"]
    vsg2["VSGSeq-individual"] --> File2["File2"]
    File2 --> Q2
    File3["File3"] --> Q2
    Q2 --> C2 & C2 & C2
    C2 -- Yes --> D2 & D2 & D2
    C2 -- No --> E2 & E2 & E2
    D2 --> F2 & F2 & F2
    E2 --> F2 & F2 & F2
    F2 --> G2 & G2 & G2
    G2 --> H2 & H2 & H2
    H2 --> I2 & I2 & I2
    I2 --> J2 & J2 & J2
    J2 --> M2 & M2 & M2
    M2 --> N2 & N2 & N2
    N2 --> P2 & P2 & P2
    P2 --> FR1["File3_results"] & FR2["File2_results"] & FR3["File3_results"]
    FR1 --> Res3["Cluster results to compare samples"]
    FR2 --> Res3
    FR3 --> Res3
    Fi1["File1"] --> M3
    vsg["vsgseq2"] --> Fi2["File2"]
    Fi2 --> M3
    Fi3["File3"] --> M3
    M3 --> C3
    C3 --> D3 & K3
    D3 -- Full --> E3
    E3 --> F3
    D3 -- Analyse --> F3
    F3 --> G3
    G3 --> H3
    H3 --> K3 & I3
    K3 --> L3
    I3 --> L3
    Fl1["File1"] --> Q1
    vsg_c["vsgseq-combined"] --> Fl2["File2"]
    Fl2 --> Q1
    Fl3["File3"] --> Q1
    Q1 --> C1
    C1 -- Yes --> D1
    C1 -- No --> E1
    D1 --> F1
    E1 --> F1
    F1 --> G1
    G1 --> H1
    H1 --> I1
    I1 --> J1
    J1 --> M1
    M1 --> N1
    N1 --> P1
     N2:::Rose
     N1:::Rose
     Q2:::Rose
     Q2:::Pine
     Q2:::Ash
     Q2:::Sky
     vsg2:::Ash
     vsg2:::Pine
     FR1:::Rose
     FR1:::Pine
     FR1:::Sky
     FR2:::Rose
     FR2:::Pine
     FR2:::Sky
     FR3:::Rose
     FR3:::Pine
     FR3:::Sky
     Res3:::Rose
     Res3:::Pine
     Res3:::Sky
     vsg:::Ash
     vsg:::Pine
     vsg_c:::Ash
     vsg_c:::Pine
    classDef Rose stroke-width:1px, stroke-dasharray:none, stroke:#FF5978, fill:#FFDFE5, color:#8E2236
    classDef Class_01 stroke:#C8E6C9
    classDef Pine stroke-width:1px, stroke-dasharray:none, stroke:#254336, fill:#27654A, color:#FFFFFF
    classDef Ash stroke-width:1px, stroke-dasharray:none, stroke:#999999, fill:#EEEEEE, color:#000000
    classDef Sky stroke-width:1px, stroke-dasharray:none, stroke:#374D7C, fill:#E2EBFF, color:#374D7C
    linkStyle 1 stroke:none,fill:none
    linkStyle 47 stroke:none,fill:none
    linkStyle 63 stroke:none,fill:none
```

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

To run the entire analysis steps, post assembly:
```
nextflow run goldrieve/vsgseq2 -r main --outdir tutorial_results --mode analyse
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