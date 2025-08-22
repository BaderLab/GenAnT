#!/bin/bash

cd $outDir

export PATH="$externalDir/my_interproscan/interproscan-5.69-101.0:$PATH"

mkdir -p interproscan ; cd interproscan

interproscan.sh -i $outDir/transcript_selection/mikado_lenient.faa -f tsv -dp
