#!/bin/bash

cd $outDir

mkdir -p interproscan ; cd interproscan

$externalDir/my_interproscan/interproscan-5.69-101.0/interproscan.sh -i $outDir/transcript_selection/mikado_lenient.faa -f tsv -dp
