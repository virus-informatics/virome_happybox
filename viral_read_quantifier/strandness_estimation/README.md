# Estimate RNA-seq strandedness from host-derived reads

## Set up a conda environment

Install [Salmon](https://github.com/COMBINE-lab/salmon) and [gffread](https://github.com/gpertea/gffread) using Conda.

```bash
conda create -n strandness_estimation -c conda-forge -c bioconda python=3.12 "salmon>=2" gffread
conda activate strandness_estimation
```

## 1. Prepare input files

The script assumes that all inputs already exist locally. The genome FASTA and annotation must be uncompressed. FASTQ files may be plain or gzip-compressed. Keep R1 and R2 in matching pair order.

For [SRR13160472](https://www.ebi.ac.uk/ena/browser/view/SRR13160472), the host is dog: the sample is labeled “MDCK infection rep1 - 6h total RNA.” If you need a reference, these commands download the dog assembly **ROS_Cfam_1.0 (GCF_014441545.1)** and its GFF3 from **NCBI Annotation Release 106**, then decompress them:

```bash
mkdir -p work/reference
reference_name=GCF_014441545.1_ROS_Cfam_1.0
reference_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/441/545/${reference_name}"

curl --fail --location --retry 3 \
  "${reference_url}/${reference_name}_genomic.fna.gz" \
  --output "work/reference/${reference_name}_genomic.fna.gz"
curl --fail --location --retry 3 \
  "${reference_url}/${reference_name}_genomic.gff.gz" \
  --output "work/reference/${reference_name}_genomic.gff.gz"

gzip -dc "work/reference/${reference_name}_genomic.fna.gz" > work/reference/host.fa
gzip -dc "work/reference/${reference_name}_genomic.gff.gz" > work/reference/host.gff
```

Both files come from the same [NCBI assembly directory](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/441/545/GCF_014441545.1_ROS_Cfam_1.0/). Their sequence identifiers must match. GFF3 must describe transcripts and their exons using `ID`/`Parent`; GTF exons must have non-empty `transcript_id` attributes. The script removes GTF records with an explicitly empty `transcript_id` before extraction.

The repository provides a random 1% subset of SRR13160472 in `../example/fastq/`; see the [dataset documentation](../example/README.md). You can also supply the full FASTQ files. The script does not download or subsample reads.

## 2. Run the workflow

Arguments are positional:

```text
bash estimate_strandness.sh SAMPLE GENOME ANNOTATION R1 R2 OUTPUT [THREADS]
```

For the prepared reference and the repository's example reads:

```bash
bash estimate_strandness.sh \
  SRR13160472 \
  work/reference/host.fa \
  work/reference/host.gff \
  ../example/fastq/SRR13160472_1.fastq.gz \
  ../example/fastq/SRR13160472_2.fastq.gz \
  work/SRR13160472 \
  8
```

Use `bash estimate_strandness.sh --help` for usage.

The script prepares a temporary annotation, extracts transcripts, builds a transcript-only index, and runs `salmon quant -l A`. The genome is symlinked into the temporary directory so its FASTA index is created there. The transcriptome and Salmon index are rebuilt on every run and deleted on exit. Set `TMPDIR` before running to choose a temporary filesystem with sufficient space.

## 3. Read the results

```bash
cat work/SRR13160472/library_type.txt
```

The output directory contains `salmon_quant/` with the original Salmon output and `library_type.txt` with the sample ID, processed read pairs, inferred library type, interpretation, compatible fragment ratio, and assigned fragments.

The interpretation follows the [Salmon library type definitions](https://salmon.readthedocs.io/en/latest/library_type.html):

| Salmon type | Interpretation |
|---|---|
| `ISR` | `fr-firststrand`: R1 is antisense to the original RNA |
| `ISF` | `fr-secondstrand`: R1 is sense to the original RNA |
| `IU` | Unstranded |
