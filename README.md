# ONT low complexity filter

## Quick start

For processing reads from fastq file:

```bash
python lc_filter.py --input example_data/reads.fastq.gz --out_pattern example_data/results
```

For extracting reads later

```bash
# just extract all good regions
seqkit subseq --bed example_data/results.bed example_data/reads.fastq.gz -j 8 | seqkit replace -p ":\." -j 40 | gzip > example_data/clean_reads.fastq.gz
# filter them from too short
filtlong --min_length 1000 --min_window_q 10 example_data/clean_reads.fastq.gz | gzip > example_data/clean_reads_filtered.fastq.gz
```

Plot regions for manual inspection

```bash
# just one read
python lc_image.py --input example_data/reads.fastq.gz --interest_id 95db81aa-9d9b-41cd-b256-21416caf0668
# all reads from file
python lc_image.py --input example_data/reads.fastq.gz --interest_file example_data/interest.txt
```

## Container creation

```bash
sudo $(which singularity) build lc_filter.sif Singularity.def
```

## Example run with singularity

```bash
singularity run lc_filter.sif --input example_data/reads.fastq.gz --out_pattern example_data/results
```
