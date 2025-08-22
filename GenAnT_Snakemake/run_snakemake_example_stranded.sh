#! /bin/bash

module load singularity
# conda activate annotation_tutorial

snakemake --configfile config_example_stranded.yaml  --jobs 750 --latency-wait 60 --cluster "qsub -cwd -V -o snakemake.stranded.output.log -e snakemake.stranded.error.log -pe smp {threads} -l h_vmem={params.memory_per_thread} {params.extra_cluster_opt} -l h_stack=32M -l h_rt={params.walltime} -P simpsonlab -b y" "$@"
