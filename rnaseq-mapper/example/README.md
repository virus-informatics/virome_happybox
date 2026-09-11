# Example Data for rnaseq-mapper

## Example dataset

| Path                                     | Description                 |
| ---------------------------------------- | --------------------------- |
| `GCF_000865725.1.sequence.fasta`         | Reference genome            |
| `SRR13160472/SRR13160472_{1,2}.fastq.gz` | Paired-end reads            |
| `sra_info.csv`                           | `--sra-info` input          |
| `strandedness_info.csv`                  | `--strandedness-info` input |

### downsampled SRR13160472

Ref: https://www.ncbi.nlm.nih.gov/sra/?term=SRR13160472

| Item          | Value                        |
| ------------- | ---------------------------- |
| LibraryLayout | `PAIRED`                     |
| strandedness  | `fr-secondstrand`            |
| Subset        | 1% of the full-size run      |
| Reads         | 152,673 per file, 51 bp each |

### GCF_000865725.1 genome fasta

Influenza A virus (A/Puerto Rico/8/1934(H1N1)). The 8 segments are manually merged into a single FASTA file.

Ref: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000865725.1/

## How to prepare these data

The FASTQ files are a random 1% subset of the full-size run, downsampled with [seqkit](https://bioinf.shenwei.me/seqkit/) v2.13.0:

```bash
seqkit sample2 --proportion 0.01 --two-pass --rand-seed 11 \
    --out-file SRR13160472_1.fastq.gz <full-size>/SRR13160472_1.fastq.gz
seqkit sample2 --proportion 0.01 --two-pass --rand-seed 11 \
    --out-file SRR13160472_2.fastq.gz <full-size>/SRR13160472_2.fastq.gz
```

- `--rand-seed 11` is shared by `_1` and `_2`, so the read pairs stay matched.
- `--two-pass` keeps the extracted read count stable for a given proportion.
- A random subset is used instead of the first N reads: this run is registered with the `alignment_sorted` schema, so its spot order may follow genomic coordinates and a head-truncated subset would be positionally biased.
