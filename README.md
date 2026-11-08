# vsgseq2

[![vsgseq2 CI](https://github.com/goldrieve/vsgseq2/workflows/vsgseq2%20CI/badge.svg)](https://github.com/goldrieve/vsgseq2/actions)

**vsgseq2** is an updated pipeline for analysing VSG-seq data, based on the original VSGSeq pipeline described in [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4514441/) and [repository](https://github.com/mugnierlab/VSGSeqPipeline).

## Installation

This pipeline is built using [Nextflow](https://www.nextflow.io/). Nextflow can be installed via Conda like this:

```bash
conda create --name nf-env bioconda::nextflow
conda activate nf-env
```

### Download test data

You will require `wget`, which can be installed using your OS-specific package manager or Conda, e.g.:

    conda install -n nf-env anaconda::wget

```bash
wget https://github.com/goldrieve/vsgseq2/raw/refs/heads/main/data/reads.tar.xz
tar -xf reads.tar.xz
cd reads
```

Run vsgseq2 on test data with Conda:

```bash
nextflow run goldrieve/vsgseq2 \
  -r main \
  --samplesheet samples.csv \
  --outdir results \
  -profile conda
```

Each module has its own Conda environment (`modules/local/<module>/environment.yml`). Use `-profile docker`, `-profile singularity` or `-profile apptainer` once module containers are available.

### Samplesheet

The samplesheet is a CSV with the columns `sample`, `r1` and `r2`. Leave `r2` empty for single-end reads. Sample names must be unique and cannot contain spaces or full stops. See `assets/samplesheet.csv`.

```csv
sample,r1,r2
early_1,early_1_1.fastq.gz,early_1_2.fastq.gz
latesingle,latesingle.fq.gz,
```

### Test profile

From a clone of the repository, run the bundled test data with:

```bash
nextflow run main.nf -profile test,conda --outdir results/test
```

## Example Output

After running the tutorial data, output files are organised into the following directories:

### 1. `VSGs/`

- `{sample}_ORF_VSGs.fasta` – Predicted VSG ORF sequences per sample
- `concatenated_vsgs.fasta` – All assembled VSGs from all samples
- `VSGome/`

  - `VSGome.fasta` – Clustered VSG database using cd-hit-est
  - `VSGome.fasta.clstr` – Summary of clustering results

### 2. `assemblies/`

- `{sample}_trinity.Trinity.fasta` – De novo transcriptome assemblies

### 3. `summary/`

- `cluster/`

  - `filtered_tpm_clusters.csv` – TPM values with cluster assignments
  - `filtered_tpm_clusters_length.csv` – As above with sequence lengths
  - `cluster_champion.csv` – Representative VSGs per cluster
  - `champion_vsgs.fasta` – Sequences of representative VSGs

- `length/`

  - `length.csv` – Sequence length data

- `multiqc_report.html` – Summary of quantification (Salmon)
- `multiqc_data/` – MultiQC data files
- `read_counts/`

  - `num_reads.csv` – Read counts per VSG
  - `total_read_counts.csv` – Total reads per sample

- `tpm/`

  - `tpm.csv` – Raw TPM values
  - `filtered_tpm.csv` – TPM values filtered by read threshold
  - `cluster_tpm.csv` – Cluster-level TPMs

- `vsgs/`

  - `vsg_count.csv` – Number of unique VSGs per sample

### 4. `trimmed/`

- `{sample}_trimmed.fq.gz` – Quality- and adapter-trimmed reads

### 5. `pipeline_info/`

- Execution report, timeline, trace, DAG and software versions

### 6. Figures

- To generate a summary figure from `cluster_tpm.csv`, run the R script `bin/plot_script.R`. It should produce this plot:

![tutorial\_figure](figures/vsg_summary.png)

## Pipeline Overview

Below is a simplified DAG of the vsgseq2 pipeline:

```mermaid
flowchart TD
    subgraph vsgseq2["."]
        C3["Trim reads"]
        D3{"Full or Analyse?"}
        E3["Trinity assembly"]
        F3["Find ORFs in contigs"]
        G3["BLAST merged ORFs vs VSG DB"]
        H3["Concatenate VSGs"]
        I3["Cluster with cd-hit-est"]
        K3["Quantify with Salmon"]
        L3["Write final results"]
        M3["Multiple files are analyzed in parallel"]
    end

    Fi1["File1"] --> M3
    Fi2["File2"] --> M3
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
```

## Running Specific Pipeline Sections

Use the `--mode` flag to control which parts of the pipeline are executed.

- `full` (default): Run the entire pipeline from raw reads
- `analyse`: Use pre-assembled transcripts to rerun downstream analyses

### Example

Run the full pipeline:

```bash
nextflow run goldrieve/vsgseq2 \
  -r main \
  --mode full \
  --samplesheet samples.csv \
  --outdir results \
  -profile conda
```

Re-use Trinity assemblies and re-run the analysis section with a new threshold:

```bash
nextflow run goldrieve/vsgseq2 \
  -r main \
  --mode analyse \
  --samplesheet samples.csv \
  --threshold 200000 \
  --outdir results/tutorial_200000 \
  -profile conda
```

## Command-line Options

Run `nextflow run goldrieve/vsgseq2 --help` for all options.

**Input/output:**

```text
--samplesheet     Samplesheet CSV (sample,r1,r2)                    [required]
--outdir          Output directory                                  [default: results/<today>]
--mode            full or analyse                                   [default: full]
--assemblies      Assemblies used in analyse mode                   [default: data/tutorial_assemblies/*_trinity.Trinity.fasta]
```

**References:**

```text
--vsg_db          VSG BLAST database                                [default: data/blastdb/concatAnTattb427.fa]
--notvsg_db       Non-VSG BLAST database                            [default: data/blastdb/vsgseq2NOTvsgs.fa]
--full_vsg_db     Additional VSGs to add to the VSGome              [default: none]
--genome          Genome used as Salmon decoys                      [default: data/genome/TriTrypDB-67_TbruceiEATRO1125_Genome.fasta]
--decoys          Salmon decoy names                                [default: data/genome/TriTrypDB-67_TbruceiEATRO1125_Genome_decoys.txt]
```

**Analysis:**

```text
--cdslength       Minimum ORF length (aa)                           [default: 300]
--partial         Allow partial ORFs                                [default: false]
--cdhit_id        cd-hit-est identity threshold (0-1)               [default: 0.94]
--cdhit_as        cd-hit-est alignment coverage (0-1)               [default: 0.94]
--threshold       Minimum reads to keep a sample in the TPM filter  [default: 100000]
```

**Resources:**

```text
--cores           Threads for Trimmomatic, Trinity and Salmon       [default: 4]
--trinitymem      Memory for Trinity (GB)                           [default: 20]
--forks           Maximum parallel tasks per process                [default: 5]
```

## Pipeline structure

vsgseq2 follows the nf-core template:

```text
main.nf                    entry point
workflows/vsgseq2.nf       main workflow
subworkflows/              initialisation, validation and completion
modules/local/<module>/    one process per module, with environment.yml and nf-test
conf/                      base resources, module options and the test profile
assets/schema_input.json   samplesheet schema
nextflow_schema.json       parameter schema
bin/                       helper scripts used by the modules
tests/                     integration test and test data
```

## Known Issues & Troubleshooting

1. **Platform Compatibility**:
   vsgseq2 has been tested on:

   - **Eddie (University of Edinburgh cluster)** running Rocky Linux 9
   - **MARS (University of Glasgow cluster) running Rocky Linux 8.10**
   - **MacOS Sequoia (Apple M1)**
   - **MacOS Sequoia (Intel)**

If you have Conda/Docker installed, it should work on most systems. If you encounter issues, please raise an issue or get in touch.

2. **Memory Usage (Trinity)**:
   Trinity is memory-hungry. If vsgseq2 stalls during assembly, try reducing resource usage:

   ```bash
   --cores 1 --trinitymem 10
   ```
