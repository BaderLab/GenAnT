#!/bin/bash

wd=$outDir/transcript_selection

cd $wd

ASSEMBLY=assembly.softmasked.fa
ASSEMBLYDIR=$outDir/assembly

UNIDB=$dataDir/uniprot_sprot

# singularity config
MIKADO_SIF=$externalDir/mikado_baderlab.sif
SINGULARITY_CACHEDIR=$wd/mikado/cache
SINGULARITY_TMPDIR=$wd/mikado/tmp
mkdir -p $SINGULARITY_CACHEDIR
mkdir -p $SINGULARITY_TMPDIR

for i in `ls -d blast/*/` ; do cat $i/blast_results.tsv >> blast/blast_results.tsv ; done

# mikado configure commads 

singularity exec --bind ${outDir},${wd},${ASSEMBLYDIR},${UNIDB} ${MIKADO_SIF} mikado serialise \
 -p 20 --start-method spawn \
 --orfs $outDir"/transcript_selection/transdecoder/transdecoder/mikado_prepared.fasta.transdecoder.bed" \
 --transcripts $outDir"/transcript_selection/mikado_prepared.fasta" \
 --tsv $outDir"/transcript_selection/blast/blast_results.tsv" \
 --json-conf $outDir"/transcript_selection/config.yaml" \
 --genome_fai $outDir"/assembly/assembly.softmasked.fa.fai" \
 -od $outDir"/transcript_selection" \
 --log $outDir"/transcript_selection/mikado_serialise.log" \
 --blast-targets $UNIDB"/uniprot_sprot_nodup.fasta" \
 --max-target-seqs 5 \
 --junctions $outDir"/transcript_selection/3-filt/portcullis_filtered.pass.junctions.bed" \
 --blast-loading-debug

singularity exec --bind ${outDir},${wd},${ASSEMBLYDIR},${UNIDB} ${MIKADO_SIF} mikado pick \
 --json-conf $outDir"/transcript_selection/config.yaml" \
 -db $outDir"/transcript_selection/mikado.db" \
 --mode lenient \
$outDir"/transcript_selection/mikado_prepared.gtf" \
 --log $outDir"/transcript_selection/mikado_pick.log"  \
 --loci-out $outDir"/transcript_selection/mikado_lenient.gff" \
 --fasta $outDir"/assembly/assembly.softmasked.fa" \
 --no-purge

ASMPATH=$outDir"/assembly/assembly.softmasked.fa"
gffread -y mikado_lenient.faa -g $ASMPATH mikado_lenient.gff



