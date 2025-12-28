# If GenAnT crashes

GenAnT is a meta-pipeline that integrates many tools and may fail for one reson or another. This markdown describes which intermediate files should be deleted before restarting the run, and why.

Here, we assume that $outDir is wherever your run of GenAnT is, either from using the scripts provided or the snakemake.

### cactus

Clearing intermediate files that will stop the job (and files that may get made even if a job failed) from cactus.

If the "jobstore" in your cactus directory exists, the job did not run to completion. Re-running the same script will fail due to a jobstore already being present. As a result, you need to clean these jobstores.

```
rm -r $outDir/cactus_alnr1/jobstore
rm -r $outDir/cactus_alnr2/jobstore
```

The hal to chain conversion may occur with an empty input file. In this instance, a chain file is formed without any acutal alignment info. These files are very small (typicall >1kb). The below script checks the filesize of the chain file and deletes it if it's under 500kb

```
find $outDir/cactus_alnr1/target-to-ref.chain -type f -size -500k -delete
find $outDir/cactus_alnr2/target-to-ref.chain -type f -size -500k -delete
```

### othrofinder

Orthofinder has (as far as I can tell) a non-optional feature where it prints the date of the orthofinder run. We sanitize this output so that if orthofinder runs correctly this shouldn't be a problem, but we have experienced instances where the snakemake gets confused. When this happens deleting the orthofinder directory and rerunning the rule does the trick.

```
rm -r $outDir/orthofinder 
```

### liftoff

liftoff makes a `db` file called `referencegff_db`. If you try to rerun liftoff (or the snakemake workflow) without deleting this file.

```
rm $outDir/liftoff/referencegff_db
```

### braker

I once saw an error where braker died while symbolic links (created by braker) still exists. This happened because our cluster went on maintenance and killed every job. Braker died when trying to re-run while these symbolic links still exist. Braker also starts from the beginning every time regardless.

```
rm -r $outDir/braker_sr
rm -r $outDir/braker_lr
```

# Empty gff files

This is something to check if GenAnT both crashes or succeeds (check post_flight_checks.md). Because there are a lot of options of different input files (i.e., number of reference species, RNA-seq data, etc.) there may be instances where an individual annotation tool fails but the rest of pipeline continues.


If the script failed before the `mikado_prepare` step: Delete empty GFF files to ensure their generation is re-run.

```
for i in $outDir/transcript_selection/*gff ; do  find $outDir/cactus_alnr2/target-to-ref.chain -type f -size -50k -delete ; echo $i ; done
```

If GenAnT failed after `mikado_prepare`, there is a manual step.

```
# Check to see which files were excludeds
ls $outDir/transcript_selection/excluded_input
```

For example I had a run where I only used one reference species for TOGA and I had both RNA-seq and ISO-seq data, so the directory looked like this:
```
braker.noRNA.gffread.gff  custom.gffread.gff  toga.r2.gffread.gff
```
Assuming `$outDir/transcript_selection/mikado_lenient.gff` was generated and the excluded gff files make sense from your input data, then you do not need to remove anything from `transcript_selection` as step 3 does not have to be rerun.

If, however, there is a gff file that you need in here (e.g., you had a second reference species for TOGA, then all of step3, some of step 4, and 5 need to be redone after generating to TOGA gene models). this is because the gene models from `toga.r2` were not included in transcript selection, and therefore all downstream analysis aside from `mirmachine` and `blast+infernal` are influenced.

If this is the case, you should:

Remove excluded gff files. This shouldn't impact the snakemake but will keep the workfow clean.

```
rm -r $outDir/transcript_selection/excluded_input/*.gff 

```

Remove the files involved in transcript selection

```
rm -r $outDir/transcript_selection/mikado* 
```

Remove the blast directory

```
rm -r $outDir/transcript_selection/blast
```

Remove transdecoder (ORF-finding) directory

```
rm -r $outDir/transcript_selection/transdecoder
```

Since the junction information comes directly from the RNA-seq/ISO-seq data, you don't actually have to remove those directories. Alternatively the vast majority of the time spent is merging the bam files, so it may be worth deleting the portcullis directories to be safe.

```
rm -r $outDir/transcript_selection/1-prep
rm -r $outDir/transcript_selection/2-junc
rm -r $outDir/transcript_selection/3-filt

```

Lastly, if this happened GenAnT could run to the end just with a set of transcripts missing (e.g., `toga.r2.gffread.gff`). If you see that an annotator in Step 2 failed but there is a final GFF file, you'll want to delete the outputs of step3 and step4.

Remove the following directories:

```
rm -r $outDir/ncRNA_analysis
rm -r $outDir/orthofinder
rm -r $outDir/gene_symbol
```





