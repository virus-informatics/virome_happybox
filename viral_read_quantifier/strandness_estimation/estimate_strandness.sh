#!/usr/bin/env bash
# Salmon 2.x workflow; index success is determined by the command exit status.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: bash estimate_strandness.sh SAMPLE GENOME ANNOTATION R1 R2 OUTPUT [THREADS]

  GENOME       Uncompressed host genome FASTA
  ANNOTATION   Uncompressed GFF3 or GTF with transcript/exon annotations
  R1, R2       Paired-end FASTQ files (plain or gzip), in matching pair order
  OUTPUT       Directory for salmon_quant/ and library_type.txt
  THREADS      Defaults to NSLOTS, or 8

Requires Salmon 2.x (salmon), gffread, python3, and awk on PATH.
Existing results are replaced after successful processing.
Temporary files are removed on exit. Set TMPDIR to choose temporary storage.
EOF
}

if [[ ${1:-} == -h || ${1:-} == --help ]]; then
  usage
  exit 0
fi
if [[ $# -lt 6 || $# -gt 7 ]]; then
  usage >&2
  exit 1
fi

sample=$1
genome_fasta=$2
annotation=$3
read1=$4
read2=$5
out_dir=$6
threads=${7:-${NSLOTS:-8}}

fail() { echo "ERROR: $*" >&2; exit 1; }
[[ $threads =~ ^[1-9][0-9]*$ ]] || fail "THREADS must be a positive integer"
[[ -n $sample && -n $out_dir ]] || fail "SAMPLE and OUTPUT must not be empty"

for command in salmon gffread python3 awk; do
  command -v "$command" >/dev/null 2>&1 || fail "Required command is unavailable: $command"
done
for input in "$genome_fasta" "$annotation" "$read1" "$read2"; do
  [[ -f $input && -s $input ]] || fail "Input file does not exist or is empty: $input"
done

# Keep the genome's .fai, transcriptome, and index away from the input files.
tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/strandness.XXXXXX")
trap 'rm -rf -- "$tmp_dir"' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
genome_absolute=$(python3 -c 'import pathlib, sys; print(pathlib.Path(sys.argv[1]).resolve())' "$genome_fasta")
ln -s "$genome_absolute" "$tmp_dir/genome.fa"

echo "Sample: $sample"
echo "Threads: $threads"
echo "Temporary directory: $tmp_dir"
echo "Output directory: $out_dir"

# NCBI GTF may contain gene records with transcript_id "". Leave GFF3
# attributes intact; stop at an optional embedded FASTA section.
echo "Preparing gene annotation"
awk -F '\t' '
/^##FASTA/ { exit }
/^#/ { print; next }
NF < 9 { next }
$9 ~ /transcript_id[[:space:]]+""/ { removed++; next }
{ print }
END { print "Removed records with empty transcript_id: " (removed + 0) > "/dev/stderr" }
' "$annotation" > "$tmp_dir/annotation.filtered.gff"

valid_exons=$(awk -F '\t' '
$3 == "exon" && ($9 ~ /transcript_id[[:space:]]+"[^"]+"/ || $9 ~ /(^|;)Parent=[^;]+/) { count++ }
END { print count + 0 }
' "$tmp_dir/annotation.filtered.gff")
[[ $valid_exons -gt 0 ]] || fail "No exons with transcript_id (GTF) or Parent (GFF3) were found"
echo "Exon records: $valid_exons"

echo "Extracting host transcripts"
gffread "$tmp_dir/annotation.filtered.gff" \
  -g "$tmp_dir/genome.fa" -w "$tmp_dir/transcripts.fa"
[[ -s $tmp_dir/transcripts.fa ]] || fail "Transcriptome FASTA was not created"
transcript_count=$(grep -c '^>' "$tmp_dir/transcripts.fa" || true)
[[ $transcript_count -gt 0 ]] || fail "No transcript sequences were extracted"
echo "Extracted transcripts: $transcript_count"

echo "Building Salmon index"
salmon index -t "$tmp_dir/transcripts.fa" -i "$tmp_dir/index" -p "$threads"
# set -e stops on index failure. Salmon 2.x uses a different index format,
# so do not require the versionInfo.json file used by Salmon 1.x.

echo "Mapping reads with automatic library-type detection"
salmon quant -i "$tmp_dir/index" -l A \
  -1 "$read1" -2 "$read2" -p "$threads" \
  -o "$tmp_dir/salmon_quant"
[[ -s $tmp_dir/salmon_quant/lib_format_counts.json ]] || fail "Salmon did not create lib_format_counts.json"

# Read the library type, falling back to meta_info.json when necessary.
# Unknown types produce an informational message; no confidence cutoff is used.
python3 - "$tmp_dir/salmon_quant" "$tmp_dir/library_type.txt" "$sample" <<'PY'
import json
import sys
from pathlib import Path

quant, report_path = Path(sys.argv[1]), Path(sys.argv[2])
counts = json.loads((quant / "lib_format_counts.json").read_text())
meta_path = quant / "aux_info/meta_info.json"
meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
expected = counts.get("expected_format")
if not expected:
    library_types = meta.get("library_types", [])
    if library_types:
        expected = library_types[0]

interpretations = {
    "ISR": "fr-firststrand (reverse-stranded; R1 antisense to the original RNA)",
    "ISF": "fr-secondstrand (forward-stranded; R1 sense to the original RNA)",
    "IU": "unstranded",
}
interpretation = interpretations.get(
    expected,
    "A less common Salmon library type was inferred; inspect lib_format_counts.json.",
)
report = (
    f"Sample: {sys.argv[3]}\n"
    f"Processed read pairs: {meta.get('num_processed', 'NA')}\n"
    f"Salmon inferred library type: {expected}\n"
    f"Interpretation: {interpretation}\n"
    f"Compatible fragment ratio: {counts.get('compatible_fragment_ratio', 'NA')}\n"
    f"Assigned fragments: {counts.get('num_assigned_fragments', 'NA')}\n"
)
report_path.write_text(report)
print(report, end="")
PY

# Replace the final results only after the analysis and report have succeeded.
mkdir -p "$out_dir"
rm -rf -- "$out_dir/salmon_quant"
mv "$tmp_dir/salmon_quant" "$out_dir/salmon_quant"
mv -f "$tmp_dir/library_type.txt" "$out_dir/library_type.txt"
echo "Salmon output: $out_dir/salmon_quant"
echo "Library-type summary: $out_dir/library_type.txt"
