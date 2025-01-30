# 4. Annotating non-coding RNA genes

Non-coding RNAs do not contain highly conserved exons and protein domains typically seen in mammalian protein-coding genes. Accordingly, identifying non-coding genes requires algorithms that do not rely on the same genomic features used in the gene-model identification algorithms described in steps 1-3 (e.g. ORF evaluation, intron-exon ratio etc.). Instead of the evaluation of ORFs to determine if the coding-gene model is functional, non-coding gene models are evaluated for their potential functionality based on whether the predicted secondary structure of that non-coding RNA matches a previously identified secondary structure.

### Identifying short non-coding structures with Infernal

Many non-coding annotations can be generated using [the RNA family (Rfam) database](https://rfam.org/), an open-access, and maintained database of non-coding RNAs. The primary tool used in non-coding gene annotation and classification is [INFERence of RNA ALignment (Infernal)](http://eddylab.org/infernal/). Briefly, Infernal builds covariance models of RNA molecules, to incorporate sequence homology and predicted RNA secondary structure in the annotation and classification on non-coding molecules in the genome. To reduce the runtime and memory requirement of this process, researchers typically pre-select sequences (seed) based on sequence homology to a non-coding database, RNA-seq alignments, and regions identified as “non-coding” in GFF post processing algorithms (e.g. Mikado).

#### Seeding by BLASTing against Rfam

To perform one round of seeding (i.e. identifying genomic regions likely to contain non-coding RNA molecules), you can first BLAST your unmasked genome FASTA file against the Rfam database of non-coding RNAs. The step-by-step process is explained [here](https://docs.rfam.org/en/latest/sequence-extraction.html), but we've also outlined our process. Start by downloading and unzipping the Rfam database:

```
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/Rfam.fa.gz
gunzip Rfam.fa.gz
```

Then filter duplicate non-coding RNAs stored within the database, as duplicates can interrupt generating the BLAST database that you'll need to build in order to BLAST against it. There are many ways of doing this; we used [seqtk](https://github.com/lh3/seqtk).

```
seqkit rmdup Rfam.fa > Rfam.rmdup.fa
```

After the duplicates are removed, you can now create a BLAST database which requires the deduplicated Rfam database as input. `-input_type fasta` specifies that the database is a FASTA file; `-dbtype nucl` indicates that the database is made of nucleotides; `-title Rfam_ncRNA` is a recognizable title for the database; `-parse_seqids` is required to keep the original sequence identifiers; `-out Rfam_ncRNA` is the output base name for the database, and is often the same as the title. This is another command that won't work if you have spaces in your working directory.

```
makeblastdb -in Rfam.rmdup.fa \
-input_type fasta \
-dbtype nucl \
-title Rfam_ncRNA \
-parse_seqids \
-out Rfam_ncRNA
```

Next, BLAST your unmasked genome FASTA file against the Rfam BLAST database. `-db Rfam_ncRNA` points to the base name of the BLAST database; `-query genome.fasta` points to your genome sequence; `-evalue 1e-6` describes the number of hits expected by chance; `-max_hsps 6` indicates a maximum of 6 alignments for any query-subject pair; `-max_target_seqs 6` indicates that a maximum of 6 aligned sequences are to be kept; `-outfmt 6` specifies the type of output from BLAST; `-out assembly.rfam.blastn` is the name of the output file. This command outputs all of the BLAST alignments found in the search that match all of the given criteria.

```
blastn -db Rfam_ncRNA \ 
-query genome.fasta \
-evalue 1e-6 -max_hsps 6 -max_target_seqs 6 -outfmt 6 \
-num_threads number_of_threads \
-out assembly.rfam.blastn
```

After this, convert the BLAST output to a BED file by extracting the chromosome, start coordinate, end coordinate, and name columns (columns 1, 7, 8,and 2 of a blastn output tsv):

```
awk -F "\t" '{print $1 "\t" $7 "\t" $8 "\t" $2}' assembly.rfam.blastn > assembly.rfam.bed
```

#### Seeding with previously identified non-coding RNA gene models

GFF files processed with Mikado will clearly indicate whether an RNA feature is non-coding. These features can be extracted from the GFF file and saved as a BED file with the same chromosome, start coodinate, end coordinate, and gene name columns. These coordinates will proceed to be added to `assembly.rfam.bed`.

To do this, first extract the non-coding RNA from Mikado:

```
grep -P "\tncRNA\t" mikado_annotation.gff > mikado_annotation.ncRNA.gff
```

Then use GFFRead to convert the GFF file to a BED file, only keeping the first four columns:

```
gffread mikado_annotation.ncRNA.gff --bed | cut -f1-4 > mikado_annotation.ncRNA.bed
```

Non-coding features can also be found from the output of LiftOff that may not have made it through Mikado's filtering steps (it doesn't matter if they did and there is overlap). Isolate non-coding features from the ouptut of LiftOff (in this case, we assume that the LiftOff output is derived from a RefSeq GFF file and can grab `gbkey=ncRNA` for all non-coding RNAs):

```
grep "gbkey=ncRNA" liftoff_annotation.gff | grep -P "RNA\t" > liftoff_annotation.ncRNA.gff
```

Again, use GFFRead to convert the GFF file to a four-column BED file:

```
gffread liftoff_annotation.ncRNA.gff --bed | cut -f1-4 > liftoff_annotation.ncRNA.bed
```

Now, isolate non-coding features from StringTie that may not have been included in the Mikado output. Since StringTie doesn't specify which features are coding and which are non-coding, we are going to assume that any features that were not included in the final output of Mikado may be non-coding. To find these excluded features, we can use the BEDTools tool, `bedtools subtract`, to remove Mikado features from the StringTie features. The `-A` flag is included to remove the entirety of the feature in the StringTie output if any of it is found in the Mikado output. Only keep transcript features rather than both transcripts and exons (specified by the grep statement following the BEDTools command).

```
bedtools subtract -A -a stringtie_annotation.gtf -b mikado_annotation.gff | grep -P "\ttranscript\t" > stringtie_annotation.ncRNA.gtf
```

Convert to a four-column BED file with GFFRead:

```
gffread stringtie_annotation.ncRNA.gtf --bed | cut -f1-4 > stringtie_annotation.ncRNA.bed
```

Finally grab a BED file of repeats from the output of Earl Grey, as these are also candidate non-coding regions for Infernal. Earl Grey already spits out a BED file that you can find in the "summaryFiles" directory (e.g. `earl_grey/species_EarlGrey/species_EarlGrey_summaryFiles/species.filteredRepeats.bed`). Isolate the first four columns of this file.

```
cut -f1-4 species.filteredRepeats.bed > earlGrey_annotation.ncRNA.bed
```

All of these BED files can now be combined.

#### Combine different seeding results

Combine candidate noncoding RNA containing genome coordinates (seeds) from each method into a master-list using concatenate:

```
cat assembly.rfam.bed \
 mikado_annotation.ncRNA.bed \
 liftoff_annotation.ncRNA.bed \
 stringtie_annotation.ncRNA.bed \
 earlGrey_annotation.ncRNA.bed > assembly_ncRNA_seed.bed
```

Sort the full BED file by coordinate:

```
bedtools sort -i assembly_ncRNA_seed.bed > assembly_ncRNA_seed.s.bed
```

At this point, you should have a BED file with four columns containing many different genomic coordinates. The first column is the contig ID, the second column is the starting base pair position, the third column is the ending base pair position, and the fourth column is the feature ID from the source it was derived. To remove some redundancy, remove rows that share the exact same contig ID, and base pair coordinates:

```
awk -F'\t' '!seen[$1 FS $2 FS $3]++' assembly_ncRNA_seed.s.bed | sponge assembly_ncRNA_seed.s.bed
```

If not done yet, index the unmasked genome FASTA file using SAMtools:

```
samtools faidx genome.fasta
```

Isolate DNA from the unmasked FASTA file that matches the genomic coordinates of the non-coding RNA seeds using `bedtools getfasta`. `-fi genome.fasta` specifies the input genome; `-bed assembly_ncRNA_seed.s.bed` are the bed coordinates that are used to specify the regions from the FASTA file to extract the sequences from; `-fo assembly_ncRNA_seed.fasta` is the name of the output FASTA file. The sequences in this FASTA file will be used by Infernal to determine ncRNA identities.

```
bedtools getfasta -fi genome.fasta -bed assembly_ncRNA_seed.s.bed -fo assembly_ncRNA_seed.fasta
```

#### Using Infernal to annotate non-coding genes

Once the probably DNA sequences that encode ncRNAs have been found (now located in `assembly_ncRNA_seed.fasta`), these can be searched against a database of covariance models (i.e. statistical models of RNA secondary structure and sequence consensus). Searching against these covariance models will help determine the identity of the non-coding genes. First, download the database of covariance models from Rfam, and unzip the file:

```
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
```

You will also need information from Rfam.clanin, which lists which models belong to the same "clan" (i.e. a group of homologous models, like LSU rRNA archaea and LSU rRNA bacteria). Download and compress this file.

```
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
cmpress Rfam.cm
```

Then, run Infernal using the cmscan function, which searches each sequence against the covariance model database. `-Z 1` indicates that the e-values are calculated as if the search space size is 1 megabase; `--cut_ga` is a flag to turn on using the GA (gathering) bit scores in the model to set inclusion thresholds, which are generally considered reliable for defining family membership; `--rfam` is a flag for using a strict filtering strategy for large databases (> 20 Gb) which accelerates the search at a potential cost to sensitivity; `--nohmmonly` specifies that the command must use the covariance models; `--tblout assembly_genome.tblout` is the output summary file of hits in tabular format; `-o assembly_genome.cmscan` is the main output file; `--verbose` indicates to include extra statistics in the main output; `--fmt 2` adds additional fields to the tabular output file, including information about overlapping hits; `--clanin Rfam.clanin` points to the clan information file; the final two positional arguments, `Rfam.cm` and `assembly_ncRNA_seed.fasta`, point to the covariance model database and FASTA file of sequences respectively. This step will take quite a while.

```
cmscan --cpu number_of_threads -Z 1 \
 --cut_ga --rfam --nohmmonly \
 --tblout assembly_genome.tblout \
 -o assembly_genome.cmscan \
 --verbose --fmt 2 \
 --clanin Rfam.clanin \
 Rfam.cm assembly_ncRNA_seed.fasta
```

Finally, the tabular output of infernal can be converted to a GFF file. The [perl script](https://raw.githubusercontent.com/nawrockie/jiffy-infernal-hmmer-scripts/master/infernal-tblout2gff.pl) to convert this output can be found in the Infernal documentation. The script can be run as follows with `--fmt2` and `--cmscan` indicating that the output of Infernal was generated with the `--fmt 2` option by cmscan. `assembly_ncRNA_seed.tblout` is the output of cmscan and the results are stored in `infernal.gff`.

```
perl infernal-tblout2gff.pl --cmscan --fmt2 assembly_genome.tblout > infernal.gff
```

### Identifying miRNAs with MirMachine

Micro RNAs (miRNAs) are not found with Infernal using the above steps, but can instead be identified using [MirMachine](https://github.com/sinanugur/MirMachine). MirMachine also relies on Infernal, but has clade-specific miRNA-specific secondary structures obtained from [MirGeneDB](https://mirgenedb.org/). To run, MirMachine needs to know the clade (`-n Mammalia`); the species name indicated by `-s`; and the softmasked genome FASTA sequence (`--genome`). We will also specify the model that MirMachine is using, which in this case is "deutero" (`-m deutero`) since mammals are within the group of deuterostome animals. Note that MirMachine is a snakemake pipeline, and if snakemake isn't installed, a cryptic error will be thrown.

```
MirMachine.py -n Mammalia -s name_of_species --genome genome.softmasked.fasta -m deutero --cpu number_of_threads
```

MirMachine outputs a GFF file with high confidence results found in `results/predictions/filtered_gff/name_of_species.PRE.gff`.

