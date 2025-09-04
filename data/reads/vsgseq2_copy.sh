#!/bin/bash

# Execute in the current working directory
#$ -cwd
#$ -N vsgseq2
#$ -V
#$ -pe sharedmem 10
#$ -l h_vmem=10G #memory per core
#$ -l h_rt=47:59:59

#set up environmental modules command
. /etc/profile.d/modules.sh

nextflow run ../../main.nf   --samplesheet samples.csv   --outdir results_full
