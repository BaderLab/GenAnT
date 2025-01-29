#!/bin/bash

wd=$outDir/transcript_selection

cd $wd

ASSEMBLY=assembly.softmasked.fa
ASSEMBLYDIR=$outDir/assembly

ASMPATH=$ASSEMBLYDIR/assembly.softmasked.fa

# singularity config
MIKADO_SIF=$externalDir/mikado_singularity/cache/mikado.2.3.5rc2.sif
SINGULARITY_CACHEDIR=$wd/mikado/cache
SINGULARITY_TMPDIR=$wd/mikado/tmp
mkdir -p $SINGULARITY_CACHEDIR
mkdir -p $SINGULARITY_TMPDIR

Rscript --vanilla $outDir/scripts/make_mikado_input.R


# mikado configure commads 


singularity exec --bind ${outDir},${wd},${ASSEMBLYDIR} ${MIKADO_SIF} mikado configure \
 --list $outDir/transcript_selection/mikado_input_sheet.txt \
 --reference $ASSEMBLYDIR/assembly.softmasked.fa \
 -y config.yaml \
 --scoring mammalian.yaml

 singularity exec --bind ${outDir},${wd},${ASSEMBLYDIR} ${MIKADO_SIF} mikado prepare \
 --json-conf config.yaml \
 --start-method spawn -p 20 \
 -od $outDir/transcript_selection
