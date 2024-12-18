# 4. Annotating non-coding RNA genes

Non-coding RNAs do not contain highly conserved exons and protein domains typically seen in mammalian protein-coding genes. Accordingly, identifying non-coding genes requires algorithms that do not rely on the same genomic features used in the gene-model identification algorithms described in steps 1-3 (e.g. ORF evaluation, intron-exon ratio etc.). Instead of the evaluation of ORFs to determine if the coding-gene model is functional, non-coding gene models are evaluated for their potential functionality based on whether the predicted secondary structure of that non-coding RNA matches a previously identified secondary structure.

Non-coding annotations can be generated using [the RNA family (Rfam) database](https://rfam.org/), an open-access, and maintained database of non-coding RNAs. The primary tool used in non-coding gene annotation and classification is [INFERence of RNA ALignment (Infernal)](http://eddylab.org/infernal/). Briefly, Infernal builds covariance models of RNA molecules, to incorporate sequence homology and predicted RNA secondary structure in the annotation and classification on non-coding molecules in the genome. To reduce the runtime and memory requirement of this process, researchers typically pre-select sequences (seed) based on sequence homology to a non-coding database, RNA-seq alignments, and regions identified as “non-coding” in GFF post processing algorithms (e.g. Mikado).

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

After the duplicates are removed, you can now create a BLAST database which requires the deduplicated Rfam database as input. `-input_type fasta` specifies that the database is a FASTA file; `-dbtype nucl` indicates that the database is made of nucleotides; `-title Rfam_ncRNA` is a recognizable title for the database; `-parse_seqids` is required to keep the original sequence identifiers; `-out Rfam_ncRNA` is the output base name for the database, and is often the same as the title.

```
makeblastdb -in Rfam_rmdup.fa \
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

GFF files processed with Mikado and RNA-seq data (and GFF files annotated from other approaches) will have “biotype” information already stored (e.g. indicating if the gene encodes ncRNA). These regions can be extracted from a GFF file and saved as a BED file with the same chromosome, start coodinate, end coordinate, and name columns. These coordinates will proceed to be added to `assembly.rfam.bed`. Here is a way to isolate these regions from a GFF file in R using the library [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html).

```
library(rtracklayer)

# Read in the GFF/GTF file
gtf <- rtracklayer::readGFF("mikado_annotation.gtf")

# Ensembl annotates biotypes. Mikado will simply return “non-coding”
# Isolate all features with the following biotypes from the GTF file
gtf_nonCoding <- gtf[ (gtf$gene_biotype %in% c("lncRNA","miRNA","rRNA","scaRNA","snoRNA",
                                  "snRNA")),]

# Only keep the transcripts that these biotypes encode
gtf_nonCoding_transcript <- gtf_nonCoding[gtf_nonCoding$type == "transcript",]

# Save new GFF file that only has non-coding transcripts
rtracklayer::export.gff3(gtf_nonCoding_transcript,"mikado_annotation_transcripts.gff3")

# Only isolate the specific columns that are needed for the BED file
gtf_nonCoding_transcript <- gtf_nonCoding_transcript[,c("seqid","start","end","transcript_id")]

# Create a BED file by exporting the object created above as a table
write.table(gtf_nonCoding_transcript, file = "mikado_annotation_noncoding.bed",quote=F,row.names = F,col.names = F,sep="\t")
```

#### Combine different seeding results

Combine candidate noncoding RNA containing genome coordinates (seeds) from each method into a master-list using concatenate:

```
cat assembly.rfam.bed mikado_annotation_noncoding.bed > assembly_ncRNA_seed.bed
```

Sort the full BED file:

```
bedtools sort -i assembly_ncRNA_seed.bed > assembly_ncRNA_seed.s.bed
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

Then, run Infernal using the cmscan function, which searches each sequence against the covariance model database. `-Z 1` indicates that the e-values are calculated as if the search space size is 1 megabase; `--cut_ga` is a flag to turn on using the GA (gathering) bit scores in the model to set inclusion thresholds, which are generally considered reliable for defining family membership; `--rfam` is a flag for using a strict filtering strategy for large databases (> 20 Gb) which accelerates the search at a potential cost to sensitivity; `--nohmmonly` specifies that the command must use the covariance models; `--tblout assembly_genome.tblout` is the output summary file of hits in tabular format; `-o assembly_genome.cmscan` is the main output file; `--verbose` indicates to include extra statistics in the main output; `--fmt 2` adds additional fields to the tabular output file, including information about overlapping hits; `--clanin Rfam.clanin` points to the clan information file; the final two positional arguments, `Rfam.cm` and `assembly_ncRNA_seed.fa`, point to the covariance model database and FASTA file of sequences respectively.

```
cmscan --cpu number_of_threads -Z 1 \
 --cut_ga --rfam --nohmmonly \
 --tblout assembly_genome.tblout \
 -o assembly_genome.cmscan \
 --verbose --fmt 2 \
 --clanin Rfam.clanin \
 Rfam.cm assembly_ncRNA_seed.fa
```

Finally, the tabular output of infernal can be converted to a GFF file. The [perl script](https://raw.githubusercontent.com/nawrockie/jiffy-infernal-hmmer-scripts/master/infernal-tblout2gff.pl) to convert this output can be found in the Infernal documentation. The script can be run as follows with `--fmt2` and `--cmscan` indicating that the output of Infernal was generated with the `--fmt 2` option by cmscan. `assembly_ncRNA_seed.tblout` is the output of cmscan and the results are stored in `assembly_ncRNA.gff`.

```
