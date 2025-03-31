# 3. Combining and filtering gene models

Completing the previous steps yields gene models from multiple homology-based annotations and transcript-assembly-based annotations. Most gene models will be identified across annotations, however some gene models will be method-specific.

#### Mikado

[Mikado](https://github.com/EI-CoreBioinformatics/mikado) is a tool designed to evaluate, combine, and filter gene models across multiple annotations in a way that mimics manual assembly curation. Mikado takes different GFF files as input, and outputs a filtered GFF file that is more accurate than any of the input annotation or evidence files on their own.

Mikado is directed by configuration and scoring files that can be customized to your annotation project. The four steps Mikado follows are:
1. Configure (Creates configuration files that guide Mikado)
2. Prepare (Prepares input GFF files for analysis)
3. Serialise (Creates database used for “Pick”)
4. Pick (Picks best transcripts for resulting annotation)

Mikado has a thorough user manual that includes a tutorial walking through how to use it: https://mikado.readthedocs.io/en/stable/Tutorial/. Before running the command to create the configuration file, it is best to organize your input GFF files that you would like to combine. The names of these files should be listed in a tab-delimited file, with one file described per row.

First column: the name of the file
Second column: short, unique identifier for the input file
Third column: True or False indicating whether or not the annotation is strand-specific
OPTIONAL COLUMNS:
Fourth column: Assignment of positive or negative weightings (e.g. if you think a particular input file is high quality, you can, say, put a 3 in that colum; a poor-quality dataset may have a -0.5)
Fifth column: True or False indicating if the annotation is a reference (important if updating an exising annotation, another function of mikado)

The easiest way to run everything is if you have all of these input files should be stored in a working directory that you are using to run Mikado, including the list of inputs. Here is an example of the tab-delimited file indicating the different input GFF files (note if items are separated by spaces instead of tabs, an error will be thrown):

```
braker.gff braker True 0 False
stringtie.gff stringtie True 0 False
liftoff.gff liftoff True 0 False
toga.gff toga True 0 False
```

#### 1. Mikado configure

To create the configuration file that runs Mikado, one must type `mikado configure` on the command line, pointing to the genome you are annotating, specifying the name of the configuration file, and pointing towards the TSV file of input files that you just created. Mikado provides a selection of scoring files you can use to cater your genome annotation to your species, which you can indicate with the `--scoring` argument - mammals will use built-in `mammalian.yaml`. An additional `--copy-scoring` flag can be used to copy a scoring file to your working directory so that you can customize it for your species. The scoring file is what Mikado will eventually use to selectively filter the input transcripts, and may need to be modified depending on what type of data or species you are working with (we’ll touch on this later). The configuration file that gets generated can also easily be modified if needed; the user can do this by either modifying their `mikado configure` command and rerunning it, or by modifying the resulting configuration file directly (e.g. using `nano conf.yaml`).

Below is an example command, with `-y` preceding the name of the configuration file, and `--reference` pointing to the soft-masked FASTA file of the genome that you're annotating.

```
mikado configure \
 --list list_of_inputs.tsv \
 --reference name_of_genome.fasta \
 -y conf.yaml \
 --scoring mammalian.yaml \
 --copy-scoring
```

#### 2. Mikado prepare

The next step is running `mikado prepare`, which requires any input GFF files you wish to combine. Since you have already created the configuration file and pointed Mikado to your list of input files, all you have to do is run `mikado prepare --json-conf conf.yaml`. This creates a GTF file containing non-redundant transcripts (`mikado_prepared.gtf`) and a corresponding FASTA file (`mikado_prepared.fasta`), as well as a log file. This step can be sped up by increasing the number of threads using the `-p` argument, and Mikado recommends adding the option `start-method spawn` when using parallelization.

```
mikado prepare \
 --json-conf conf.yaml
 --start-method spawn \
 -p number_of_threads
```

Our scripts perform these analyses with `12_mikado_configure_and_prepare.sh`. First, Our pipeline copies every `gff` file into `transcript_selection` once it is generated in "2-generating models".

This script scans through the `transcript_selection` directory for non-empty gff files. It then runs our `scripts/make_mikado_input.R` script to generate the `mikado_input_sheet.txt` described above before running `mikado configure` and `mikado prepare`

```
tutorialDir=/path-to-GAT/GenomeAnnotationTutorial
externalDir=$tutorialDir/external
outDir=/path-to-output/
customRef="none" # if you're building on a reference annotation
liftoffRef="none" # if the liftoff annotation should consider a reference genome (e.g., you are doing an assembly/annotation upgrade for a species with an existing assembly/annotation)

scripts/12_mikado_configure_and_prepare.sh \
$outDir \
$externalDir \
$customRef \
$liftoffRef

```

### Building transcript feature table

Before running `mikado serialise`, additional work should be done to provide Mikado with more information about the transcripts now stored in the `mikado_prepared` files. The steps are as follows, and are especially important when working with RNA-seq-derived gene models:
1. Validate splice junctions with [Portcullis](https://github.com/EI-CoreBioinformatics/portcullis)
2. Determine sequence similarity with [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
3. Identify open reading frames (ORFs) with [TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki)

#### Splice junctions

Splice junctions are crucial to defining CDS regions and intron-exon boundaries. In principle, a splice junction from RNA-seq is defined as a high quality read mapping to two parts of the same gene model. Extra filters (e.g., a canoncial splice site in the multimapping read) improves the confidence of a splice junction. Junctions and lncRNAs are two of the most strongly impacted features by the presence of RNAseq/ISOseq.

##### Portcullis (RNA-seq)

Validating splice junctions can be done with [Portcullis](https://github.com/EI-CoreBioinformatics/portcullis), which filters out false positive intron/exon boundaries which are often found in the outputs of RNA-seq alignment tools. Portcullis can be run on a merged BAM file of all of your aligned RNA-seq reads. So if you had aligned RNA-seq datasets with RNA-seq alignment tools, you would have ended up with a BAM file for each alignment performed. All of these BAM files must first be merged with [SAMtools](https://www.htslib.org/). Three are given as an example, but any number of BAM files can be merged.

```
samtools merge \
 -@ number_of_threads \
 merged_bams.bam \
 bam_file_1.bam \
 bam_file_2.bam \
 bam_file_3.bam
```

Portcullis can now be run on the output, `merged_bams.bam`; this tool analyses all of the splice junctions in the BAM file and filters out the junctions that are not likely to be genuine. Portcullis has three steps: `prep`, which prepares the data for junction analysis; `junc`, which calculates junction metrics; and `filt`, which separates valid and invalid splice junctions. The easiest way to run the tool, however, is to run `portcullis full` which combines all three steps and produces a BED file of junctions that can be used as input for `mikado serialise` (`portcullis.pass.junctions.bed`).

```
portcullis full \
 -t number_of_threads \
 name_of_genome.fasta \
 merged_bams.bam
```

##### regtools (ISO-seq)

Pulling splice junctions from ISO-seq data is very valuable, as a transcript contains mouse (if not all) junctions for a particular isofom. The nature of these data also makes junctions easier to find. We use `regtools` to pull splice junctions.

ISO-seq bam files are merged in the same was as RNA-seq

```
samtools merge \
 -@ number_of_threads \
 merged.lr.bam \
 bam_file_1.bam \
 bam_file_2.bam \
 bam_file_3.bam
```

Junctions are computed with

```
regtools junctions extract -s XS -o ./lr_junctions.bed merged.lr.bam
```

Assuming you have short and long read RNA-seq data, concatenate the junctions from regtools and Porticullus to get the final set.
```
cat 3-filt/portcullis_filtered.pass.junctions.bed lr_junctions.bed > junctions.final.us.bed
```

Using our scripts, junctions are computed with `scripts/get_junctions.sh`


Assuming the directory structure built in our pipeline/snakemake, getting junctions is performed with:
```
outDir=/path-to-output/
scripts/get_junctions.sh $outDir
```



#### BLAST+

[BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi) can now be used to identify sequence similarity to known proteins. Different protein databases exist against which the predicted transcript sequences output by `mikado prepare` can be compared; we used the high-quality curated protein database, SwissProt. This database can be downloaded from [uniprot.org](https://www.uniprot.org/).

`wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz`

Once the database is downloaded, you can use BLAST to index it which is required to perform a BLAST+ search. `-dbtype` indicates that this is a protein database. `uniprot_sprot` is the base name of the FASTA file and the resulting BLAST database.

```
makeblastdb \
 -in uniprot_sprot.fasta \
 -dbtype prot \
 -out uniprot_sprot
```

BLAST the transcript sequences from `mikado prepare` against the SwissProt database using `blastx`, which is the command required to compare translated nucleotide sequences to protein sequences. BLAST needs to be run requesting the following output format: `-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop"`. This creates a TSV results file that Mikado uses to filter or score transcripts based on their homology to existing protein-coding sequences. `-max_target_seqs` indicates to keep a maximum of this number of hits (we used 5, as seen on the Mikado tutorial). `-query` is the transcript file BLASTed against the protein database. `-outfmt` specifies the format required by the next step of Mikado. `-db` is the SwissProt database. `-evalue` is a minimum measure of significance to consider a protein sequence in the SwissProt database a hit against the query.

```
blastx \
 -max_target_seqs num_of_seqs \
 -num_threads number_of_threads \
 -query mikado_prepared.fasta \
 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop" \
 -db uniprot_sprot \
 -evalue 0.000001 \
 -out blast_results.tsv
```

The output is a TSV file of BLAST+ results called `blast_results.tsv`. Note that BLAST+ takes a very long time (maybe a day or so), also depending on the number of sequences output by `mikado prepare`. This time will increase significantly if a larger protein database is used. [Diamond](https://github.com/bbuchfink/diamond) is a much faster alternative than BLAST+, but it finds fewer hits even in ultra-sensitive mode.

BLAST has considerable runtime, and it scales with the number of predicted transcripts. However, each BLAST query is independent of each transcript, meaning we can split the mikado_prepared.fasta file and parallelize this process. Our tutorial splits the fasta file into 100 segments and runs them in parallel.

`bash scripts/splitfa.sh 100 mikado_prepared.fasta` splits your fasta into 100 segments and puts them in the `blast` directory.


`scripts/sub_mikado_blast.sh` is an example script of how to run the 100 blast jobs. 

E.g. (on our SGE).
```
for i in *fasta ; do qsub -N $i -P simpsonlab -cwd -V -v I=$i scripts/run_looped_mikado_blast.sh ; done
```


#### TransDecoder

Open reading frames (ORFs) can be determined with [TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki), a tool that scans transcripts for potential coding regions. A single transcript can produce multiple ORFs, and the user can set a minimum amino acid count using the `-m` flag. The default setting of TransDecoder is `-m 100` indicating that the minimum ORF detected is 100 amino acids long. Lowering this parameter increasing the number of false positive ORFs, but allows for shorter ORFs to be detected (we chose to set `-m 30`; Mikado also has an internal minimum ORF length of 50bp). The input to TransDecoder is the FASTA file generated by `mikado prepare`; `TransDecoder.LongOrfs` identifies the ORFs.

```
TransDecoder.LongOrfs \
 -t mikado_prepared.fasta \
 -m number_of_amino_acids
```

Next, run `TransDecoder.Predict` to predict the likely coding regions. The `--retain_long_orfs_length` argument retains all ORFs longer than this number of amino acids (again, we used 30).

```
TransDecoder.Predict \
 -t mikado_prepared.fasta \
 --retain_long_orfs_length number_of_bases
```

If there is an error in the second step, it may be fixed by adding the flag `--no_refine_starts`. This prevents the identification of potential start codons for 5' partial ORFs using a position weight matrix; this process may fail if there are not enough sequences to model the start site. Further, TransDecoder may fail if there are spaces in any of the file paths you are using.

TransDecoder outputs valid ORFs in a BED file (e.g. `mikado_prepared.fasta.transdecoder.bed`) that can now be used for Mikado serialise.

Using our pipeline, TransDecoder is run with `scripts/run_transdecoder.sh`

```
outDir=/path-to-output/
bash scripts/run_transdecoder.sh $outDir
```

#### 3. Mikado serialise

Now that you have found junctions (if possible), BLAST+, and TransDecoder, `mikado serialze` combines this information into an SQL database that generates the feature table. Some of the different files and parameters that should be provided to Mikado include:
1. The FASTA output of `mikado prepare` (`mikado_prepared.fasta)
2. The configuration file (`conf.yaml`)
3. The genome index file (`name_of_genome.fai`)
4. The output of Portcullis (`portcullis.pass.junctions.bed`)
5. The output of BLAST (`blast_results.tsv`)
6. The FASTA file used for the BLAST search (`uniprot_sprot.fasta`)
7. The maximum number of discrete hits that can be assigned to a single sequence (we set this to 5)
8. The output of TransDecoder (`mikado_prepared.fasta.transdecoder.bed`)
9. The name of the log file to be produced by `mikado serliase` (`mikado_serialise.log`)

Here is an example of a Mikado serialise command:

```
mikado serialise \
 -p number_of_threads --start-method spawn \
 --transcripts mikado_prepared.fasta \
 --json-conf conf.yaml \
 --genome_fai name_of_genome.fai \
 --junctions portcullis.pass.junctions.bed \
 --tsv blast_results.tsv \
 --blast-targets uniprot_sprot.fasta \
 --max-target-seqs number_of_targets \
 --orfs mikado_prepared.fasta.transdecoder.bed \
 --log mikado_serialise.log
```

This creates a database called `mikado.db`. Note that if one of your input sources changes and you want to rerun `mikado serliase`, you have to manually delete `mikado.db` or else an error will be thrown.

#### 4. Mikado pick

The final step of the Mikado pipeline, `mikado pick`, takes this file as input and selects what it determines to be the best gene models. In order to perform this `pick` command, Mikado relies on a scoring file (e.g. `mammalian.yaml`) that guides the algorithm on what parameters create the best gene models. For instance, the best transcripts may long sequences with more than two exons, and introns less than 2000bp long (more parameters are considered, that was just an example). These parameters are described in the scoring file which may be customized by the user. Instructions on how to do this can be found in the [Mikado guidelines](https://mikado.readthedocs.io/en/stable/Tutorial/Scoring_tutorial/#configure-scoring-tutorial). Mikado has a flag `--no-purge` which can be used to prevent Mikado from throwing out gene models that fail specific requirements in the scoring and configuration files, but where there is no competing gene model. This dramatically increases the number of gene models in the final annotation, and we have found that it yields higher BUSCO scores (at the risk of including more false positives).

This is also where the user has to decide how they want to treat the chimeras in their gene models. Mikado has five different options that can be chosen by the user which range in stringency: nosplit, stringent, lenient, permissive, split. “nosplit” never splits any gene models, whereas “split” always splits multi-ORF transcripts. The other options land somewhere in between and really on homology results from BLAST hits (e.g. only splitting if consecutive ORFs have BLAST hits against two different targets).

```
mikado pick \
 --json-conf conf.yaml \
 -db mikado.db \
 --mode lenient \
 mikado_prepared.gtf \
 --scoring mammalian.yaml
 --loci-out name_of_final_annotation.gff \
 --log mikado_pick.log \
 --no-purge
```

The output of `mikado pick` is a GFF file containing the gene models selected based on the parameters in the scoring file and the information in `mikado_serialise.db`. At this point, the gene models can be analysed and visualised (if desired) for quality purposes. `mikado pick` is fairly quick to run, so it may be a good idea to run it a few times using different stringency levels on when to split chimeras, to see which setting results in the most expected gene model statistics (e.g. the highest BUSCO scores).

In our scripts, serialize and pick are run using th `scripts/45_mikado_serialize_pick.sh` command.

```
tutorialDir=/path-to-GAT/GenomeAnnotationTutorial
externalDir=$tutorialDir/external
datalDir=$tutorialDir/data
outDir=/path-to-output/

scripts/45_mikado_serialize_pick.sh {config[outDir]} {config[dataDir]} {config[externalDir]}
```


#### Mikado and associated tools: installing/running/troubleshooting

- Mikado has been challenging to install as it has a lot of dependencies. Therefore, we created a Docker image that can be run as follows: `docker run -v "$(pwd)":/tmp risserlin/mikado:ubuntu22_mikado2.3.2 mikado --help`
- Mikado has many steps but should not take more than two days to run; the slowest steps are BLAST+ and TransDecoder, followed by `mikado serialise`
- TransDecoder and Portcullis can be installed with Conda
- TransDecoder and creating a blast database do not work with spaces in the file paths
- `mikado_prepared.fasta.fai` and `mikado.db` need to be manually deleted if rerunning the whole Mikado pipeline in the same directory as the files will not be overwritten and confusing errors will be thrown
- Make sure that the input list of samples is a TSV separated file; spaces separating each column will throw an error
