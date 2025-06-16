#!/bin/bash

outDir=$1
externalDir=$2
customRef=$3
liftoffRef=$4
snakeDir=$5
mikadoScore=$6 

wd=$outDir/transcript_selection

cd $wd

mkdir -p excluded_input 

echo "scanning potential gff files to be used in transcript selection".

for i in *.gffread.gff ; do 
	if [[ $(wc -l < $i) -gt 5 ]] ; then

	echo "$i is included in transcript selection."

	else

	echo "$i is excluded from transcript selection."
	mv $i excluded_input

	fi
done

ASSEMBLY=assembly.softmasked.fa
ASSEMBLYDIR=$outDir/assembly

ASMPATH=$ASSEMBLYDIR/assembly.softmasked.fa

# singularity config
MIKADO_SIF=$externalDir/singularity_images/mikado_gat.sif
SINGULARITY_CACHEDIR=$wd/mikado/cache
SINGULARITY_TMPDIR=$wd/mikado/tmp
mkdir -p $SINGULARITY_CACHEDIR
mkdir -p $SINGULARITY_TMPDIR

Rscript --vanilla $snakeDir/scripts/make_mikado_input.R -c $customRef -l $liftoffRef


# mikado configure commads 


singularity exec --bind ${outDir},${wd},${ASSEMBLYDIR} ${MIKADO_SIF} mikado configure \
 --list $outDir/transcript_selection/mikado_input_sheet.txt \
 --reference $ASSEMBLYDIR/assembly.softmasked.fa \
 -y config.yaml \
 --scoring $mikadoScore

 singularity exec --bind ${outDir},${wd},${ASSEMBLYDIR} ${MIKADO_SIF} mikado prepare \
 --json-conf config.yaml \
 --start-method spawn -p 20 \
 -od $outDir/transcript_selection

samtools faidx $outDir/transcript_selection/mikado_prepared.fasta
