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
-evalue 1e-2 -max_hsps 6 -outfmt 6 \
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

Finally grab a BED file of repeats from the output of Earl Grey, as these are also candidate non-coding regions for Infernal. Earl Grey already spits out a BED file that you can find in the "summaryFiles" directory (e.g. `earl_grey/species_EarlGrey/species_EarlGrey_summaryFiles/species.filteredRepeats.bed`). To isolate only the non-coding RNAs found by Earl Grey and not the repeats, we can do a `grep` search for "RNA", using yhr `-i` flag to ignore case restrictions in case any parts of "RNA" are present in lowercase. We'll then pipe this command into `cut` to isolate the first four columns of this file.

```
grep -i "RNA" species.filteredRepeats.bed | cut -f1-4 > earlGrey_annotation.ncRNA.bed
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

At this point, you should have a BED file with four columns containing many different genomic coordinates. The first column is the contig ID, the second column is the starting base pair position, the third column is the ending base pair position, and the fourth column is the feature ID from the source it was derived. To remove some redundancy, we can use `bedtools merge` to dissolve overlapping genomic locations, since Infernal only needs the general locations for seeding and doesn't rely on the exact coordinates for anything.

```
bedtools merge -i assembly_ncRNA_seed.s.bed | sponge assembly_ncRNA_seed.s.bed
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

Looking at `infernal.gff`, you will see that the type of ncRNA that Infernal has identified in column three. tRNAs are easily recognizable, but all other ncRNAs are labeled with a more cryptic RFam identifier, like "SSU_rRNA_bacteria", "U6", etc. In order to better interpret these results, we will want to reformat the GFF so that the feature types are more recognizable. We can do that by first downloading a table that contains the different RFam identifiers and their respective feature types.

```
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/family.txt.gz
gunzip family.txt.gz
```

The easiest way to use this table is to read it into R, in addition to the GFF file from Infernal, and use pattern matching to assign the different feature types to the ncRNAs output by Infernal. We have provided an R notebook that performs these modifications called `rfamConversion.Rmd` which can be found in this folder. The notebook outputs a reformatted GFF called `infernal.types.gff`.

### Identifying miRNAs with MirMachine

Micro RNAs (miRNAs) are not found with Infernal using the above steps, but can instead be identified using [MirMachine](https://github.com/sinanugur/MirMachine). MirMachine also relies on Infernal, but has clade-specific miRNA-specific secondary structures obtained from [MirGeneDB](https://mirgenedb.org/). To run, MirMachine needs to know the clade (`-n Mammalia`); the species name indicated by `-s`; and the softmasked genome FASTA sequence (`--genome`). We will also specify the model that MirMachine is using, which in this case is "deutero" (`-m deutero`) since mammals are within the group of deuterostome animals. Note that MirMachine is a snakemake pipeline, and if snakemake isn't installed, a cryptic error will be thrown.

```
MirMachine.py -n Mammalia -s name_of_species --genome genome.softmasked.fasta -m deutero --cpu number_of_threads
```

MirMachine outputs a GFF file with high confidence results found in `results/predictions/filtered_gff/name_of_species.PRE.gff`.

The GFF file that MirMachine is a little funky, as it has `gene_id=` instead of `ID=` which is expected for GFF3 files. We can use sed to replace `gene_id` with `ID`.

```
cd results/predictions/filtered_gff/
sed 's/gene_id/ID/g' name_of_species.PRE.gff > name_of_species.PRE.id.gff
```

Now to clearly have these features labeled as ncRNAs when they are combined into the full GFF with protein-coding gene models, let's reformat the file with GFFRead and then add `gbkey=ncRNA` onto the end of each line. We can actually use this value to replace everything including and after `sequence_with_30nt` which is currently at the end of each line and is quite verbose. We can just use `sed -i` to make the edits in place.

```
sed -i 's/sequence_with_30nt.*/gbkey=ncRNA/g' name_of_species.PRE.id.gff
```

### Combining ncRNA gene models

At this point, we now have three different sources of ncRNA from Mikado, MirMachine, and Infernal. Many of these gene models will likely be overlapping, so we will have to determine which we would like to keep. We can do that by prioritizing certain gene models over others. Let's deal with lncRNAs first. Infernal found a bunch of lncRNAs which we labeled as such. Most of Mikado's ncRNA gene models are likely lncRNAs; Mikado may also have short mRNA transcripts that have been mistakenly assigned a coding region but are actually lncRNAs. One tricky thing is that Infernal doesn't find whole lncRNA transcripts, just exons, whereas Mikado can find multiexonic transcripts. So Infernal will have labeled multiple exons of a lncRNA separately without joining them together whereas Mikado may have found the whole transcript. We therefore want to add the lncRNAs from Infernal to Mikado, and if there is any overlap in the location of the gene models, we want to keep Mikado's gene model but Infernal's ncRNA label.

First, let's isolate the lncRNA features from the Infernal output using `grep`. The `-P` flag indicates that we're using the perl language to specify that we want to use the `\t` symbol to represent tabs. We name the output `infernal.types.lncRNA.gff`.

```
grep -P '\tlncRNA\t' infernal.types.gff > infernal.types.lncRNA.gff
```

Then use `bedtools intersect` to find any lncRNAs from Infernal that don't have any overlap at all with Mikado's gene models (referred to as `mikado.gff`. `-v` indicates that we're looking for models with no overlap between the GFFs specified by `-a` and `-b`, specifically keeping the non-overlapping models from `-a`.

```
bedtools intersect -v -a infernal.types.lncRNA.gff -b mikado.gff > infernal.lncRNA.notInMikado.gff
```

If any gene models have populated `infernal.lncRNA.notInMikado.gff`, append them to `mikado.gff`.

```
cat mikado.gff infernal.lncRNA.notInMikado.gff > mikado.infernal.gff
```

Now we want to label anything in `mikado.infernal.gff` with information from Infernal. We can do this using `bedtools intersect` again with `-a` pointing at the Mikado GFF file (with the Infernal-only lncRNAs added) and `-b` pointing at all of the Infernal lncRNAs. The `-wo` flag indicates that we want to write all of the Mikado features that overlap with the Infernal features, but keeping the information from both files. Therefore all possible lncRNA positions identified by Infernal will be present in this file, but associated with the transcript coordinates identified by Mikado. We give this a text file extension because it no longer follows GFF format.

```
bedtools intersect -a mikado.infernal.gff -b infernal.types.lncRNA.gff -wo > mikado.infernal.lncRNALabeled.txt
```

Let's cut the first nine columns to extract the Mikado information, and the 10th to 18th columns to extract the Infernal information. Giving each file nine columns will return proper GFF format. We don't need to worry about the number of base pair overlaps in column 19.

```
cut -f1-9 mikado.infernal.lncRNALabeled.txt > mikado.infernal.lncRNALabeled.mikadoInfo.gff
cut -f10-18 mikado.infernal.lncRNALabeled.txt > mikado.infernal.lncRNALabeled.infernalInfo.gff
```

These two separate files can then be read into R via an R notebook that we provide called `RenameLncRNAs.Rmd` that (1) renames Mikado's gene models with Infernal lncRNA models found in the same location and (2) reformats Infernal gene models that did not overlap with Mikado's. The R notebook outputs a file called `mikado.infernal.lncRNALabeled.polished.gff` which is still a subset of lncRNAs that will need to be reintegrated with the rest of Mikado's gene models.

We can integrate these newly formatted lncRNAs with the rest of Mikado's gene models by first subtracting the lncRNA features from the original Mikado features. Specifically, we'll be removing gene models from the file we created earlier called `mikado.infernal.gff`. This will provide us with a file that only has Mikado features that don't share any overlap with those in `mikado.infernal.lncRNALabeled.polished.gff`. We use the `-A` flag which indicates that an entire feature is removed if there is any overlap to prevent features from being broken up into smaller pieces (which shouldn't matter because we only extracted whole features earlier, but it's a good precaution). We'll call this file `mikado.noLnc.gff`.

```
bedtools subtract -A -a mikado.infernal.gff -b mikado.infernal.lncRNALabeled.polished.gff > mikado.noLnc.gff
```

Now we'll add those nicely-formatted lncRNA features from `mikado.infernal.lncRNALabeled.polished.gff` back into `mikado.noLnc.gff` using the `cat` command, which will just append the features from the former file onto the latter. We won't worry about sorting the GFF file by coordinate until the short non-coding RNAs are added, as well.

```
cat mikado.noLnc.gff mikado.infernal.lncRNALabeled.polished.gff > mikado.lncLabeled.gff
```

At this point, we now want to add the short non-coding RNAs to the gene models found in `mikado.lncLabeled.gff`. There will be two key steps here: (1) Removing any non-coding RNAs that exist within an exon found by Mikado (because it means that the non-coding RNA probably isn't real), and (2) creating gene and exon features for the remaining non-coding RNAs.

Let's start by isolating all of the exons from `mikado.lncLabeled.gff` using a simple grep command. `-P` means we're using Perl to recognize tabs on either side of "exon" to isolate rows where "exon" is in the feature type column, column two. We'll need this to check for overlap later with BEDTools.

```
grep -P "\texon\t" mikado.lncLabeled.gff > mikado.lncLabeled.exons.gff
```

Then we need to remove the lncRNAs from the Infernal results since we already dealt with those. We can use grep for that, as well. The `-v` flag means we're excluding whatever we're searching for, and `-P` allows us to specify tabs on either side of "lncRNA".

```
grep -v -P "\tlncRNA\t" infernal.types.gff > infernal.types.noLncRNA.gff
```

Now we'll want to combine all of the small non-coding RNAs from Infernal and MirMachine. Because both tools find microRNAs, there may be overlapping gene models. We're going to prioritize those from MirMachine, as MirMachine has clade-specific models that can be used allowing it to be more specific with its microRNA detection. Therefore, we're going to use BEDTools similarly to above. Let's start by subtracting any MirMachine microRNA models from the Infernal short non-coding RNAs. `-A` indicates that we're subtracting entire features with any overlap, and the file indicated by `-b` is getting subtracted from the file indicated by `-a`.

```
bedtools subtract -A -a infernal.types.noLncRNA.gff -b name_of_species.PRE.id.gff > infernal.noLnc.noMir.gff
```

Now let's concatenate the file we just made with the MirMachine GFF to collect all short non-coding RNAs.

```
cat name_of_species.PRE.id.gff infernal.noLnc.noMir.gff > short_ncRNAs.gff
```

Now we will want to remove any gene models from `short_ncRNAs.gff` that overlap with `mikado.lncLabeled.exons.gff`. We can use `bedtools subtract` again to remove the gene models, using `-A` to indicate that we are removing entire features if there is any overlap.

```
bedtools subtract -A -a short_ncRNAs.gff -b mikado.lncLabeled.exons.gff > short_ncRNAs.noOverlap.gff
```

Also note that Infernal and MirMachine both have recorded "e-values", a statistic indicating the confidence of the gene model. However, MirMachine's e-value is noted as "E-value" whereas Infernal's is "evalue" in the GFF file. We can use `sed -i` to modify the file in-place, replacing "E-value" with "evalue" for cleanness and consistency.

```
sed -i 's/E-value/evalue/g' short_ncRNAs.noOverlap.gff
```

Before being added to the Mikado + lncRNA gene models, these short ncRNAs are missing gene and exon features. We can add these features via an R notebook called
