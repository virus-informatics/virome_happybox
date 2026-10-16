# Example Data for viral_read_quantifier

## Example dataset

| Path                                     | Description                 |
| ---------------------------------------- | --------------------------- |
| `GCF_000865725.1.sequence.fasta`         | Reference genome            |
| `fastq/SRR13160470_{1,2}.fastq.gz`       | Paired-end reads for SRR13160470 |
| `fastq/SRR13160472_{1,2}.fastq.gz`       | Paired-end reads for SRR13160472 |
| `sra_info.csv`                           | `--sra-info` input          |
| `strandedness_info.csv`                  | `--strandedness-info` input |

### Downsampled sample data

The example FASTQ files are random 1% subsets of these paired-end runs:

| Run                                                               | LibraryLayout | strandedness      | Reads per FASTQ file | Read length |
| ----------------------------------------------------------------- | ------------- | ----------------- | -------------------: | ----------: |
| [SRR13160470](https://www.ncbi.nlm.nih.gov/sra/?term=SRR13160470) | `PAIRED`      | `fr-secondstrand` |              146,915 |       51 bp |
| [SRR13160472](https://www.ncbi.nlm.nih.gov/sra/?term=SRR13160472) | `PAIRED`      | `fr-secondstrand` |              152,673 |       51 bp |

### GCF_000865725.1 genome fasta

Influenza A virus (A/Puerto Rico/8/1934(H1N1)). The 8 segments are manually merged into a single FASTA file.

Ref: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000865725.1/

## How to prepare these data

The FASTQ files are random 1% subsets of the full-size runs, downsampled with [seqkit](https://bioinf.shenwei.me/seqkit/) v2.13.0:

```bash
for run in SRR13160470 SRR13160472; do
    seqkit sample2 --proportion 0.01 --two-pass --rand-seed 11 \
        --out-file "fastq/${run}_1.fastq.gz" "<full-size>/${run}_1.fastq.gz"
    seqkit sample2 --proportion 0.01 --two-pass --rand-seed 11 \
        --out-file "fastq/${run}_2.fastq.gz" "<full-size>/${run}_2.fastq.gz"
done
```

- `--rand-seed 11` is shared by `_1` and `_2`, so the read pairs stay matched.
- `--two-pass` keeps the extracted read count stable for a given proportion.
- Random subsets are used instead of the first N reads: these runs are registered with the `alignment_sorted` schema, so their spot order may follow genomic coordinates and a head-truncated subset would be positionally biased.
