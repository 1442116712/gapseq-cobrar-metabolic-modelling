#!/bin/bash
set -e

TARGET=${1:?Usage: $0 <SGB_ID>}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/pipeline.config"

WORKDIR="$PIPELINE_BASE/${TARGET}_rep"
META="$PIPELINE_BASE/genomes-all_metadata.tsv"

mkdir -p "$WORKDIR"
cd "$WORKDIR"

[ -f "$META" ] || wget -q -O "$META" \
    ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/cow-rumen/v1.0.1/genomes-all_metadata.tsv

N=$(awk -F'\t' -v t="$TARGET" '$14 == t' "$META" | wc -l)
echo "Cluster $TARGET has $N members"

mkdir -p gff fna

awk -F'\t' -v t="$TARGET" '$14 == t {print $1"\t"$20}' "$META" | \
while IFS=$'\t' read -r id url; do
    out="gff/${id}.gff.gz"
    [ -s "$out" ] && continue
    wget -q -O "$out" "$url" || echo "FAIL gff $id"
done

for gff in gff/*.gff.gz; do
    id=$(basename "$gff" .gff.gz)
    [ -s "fna/${id}.fna" ] && continue
    zcat "$gff" | awk '/^##FASTA$/ {found=1; next} found' > "fna/${id}.fna"
done

for fna in fna/*.fna; do
    id=$(basename "$fna" .fna)
    [ -s "${id}.faa" ] && continue
    prodigal -i "$fna" -a "${id}.faa" -p single -q -o /dev/null
done

echo "GFF non-empty: $(find gff -name '*.gz' -size +0 | wc -l)"
echo "FNA non-empty: $(find fna -name '*.fna' -size +0 | wc -l)"
echo "FAA non-empty: $(find . -maxdepth 1 -name '*.faa' -size +0 | wc -l)"
