#! /bin/bash

module load singularity
# conda activate annotation_tutorial

snakemake --jobs 750 --latency-wait 60 --cluster "qsub -cwd -V -o snakemake.output.log -e snakemake.error.log -pe smp {threads} -l h_vmem={params.memory_per_thread} {params.extra_cluster_opt} -l h_stack=32M -l h_rt={params.walltime} -P simpsonlab -b y" "$@"
