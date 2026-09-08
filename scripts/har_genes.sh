#!/usr/bin/env bash
set -euo pipefail

# ======================================================================
# Get gene IDs for genes falling into Human Accelerated Regions (HARs)
# Default genome: hg38
#
# Outputs:
#  1) <outprefix>.genes_overlapping_HARs.bed
#       chrom  start  end  gene_id  gene_symbol  strand
#  2) <outprefix>.ensembl_gene_ids.txt (if annotation=gencode)
#  3) <outprefix>.nearest_TSS.tsv (if --mode nearest or --mode both)
#
# Usage:
#   bash get_har_genes.sh [--mode overlap|nearest|both]
#                         [--annotation gencode|refseq]
#                         [--outprefix OUTPUT_PREFIX]
#
# Requirements: bash, awk, bedtools, wget (or curl), gunzip
# ======================================================================

MODE="overlap"           # overlap | nearest | both
ANNOT="gencode"          # gencode | refseq
OUTPREFIX="har_genes"

# HARs source (hg38): from GEO GSE180714 supplementary files
HAR_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180714/suppl/GSE180714%5FHARs.bed.gz"

# GENCODE v44 (hg38) GTF
GENCODE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"

# RefSeq (hg38) table from UCSC
REFSEQ_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz"

# -------------------- arg parsing --------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode)
      MODE="${2}"; shift 2;;
    --annotation)
      ANNOT="${2}"; shift 2;;
    --outprefix)
      OUTPREFIX="${2}"; shift 2;;
    -h|--help)
      sed -n '1,70p' "$0"; exit 0;;
    *)
      echo "Unknown arg: $1" >&2; exit 1;;
  esac
done

if [[ "$MODE" != "overlap" && "$MODE" != "nearest" && "$MODE" != "both" ]]; then
  echo "ERROR: --mode must be overlap|nearest|both" >&2; exit 1
fi

if [[ "$ANNOT" != "gencode" && "$ANNOT" != "refseq" ]]; then
  echo "ERROR: --annotation must be gencode|refseq" >&2; exit 1
fi

# -------------------- dependency checks --------------------
need() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH"; exit 1; }
}
need awk
need gunzip
need bedtools

# If wget not present, try curl
DL="wget -O"
if ! command -v wget >/dev/null 2>&1; then
  if command -v curl >/dev/null 2>&1; then
    DL="curl -L -o"
  else
    echo "ERROR: need wget or curl" >&2; exit 1
  fi
fi

# -------------------- download HARs --------------------
HAR_GZ="GSE180714_HARs.bed.gz"
HAR_BED="GSE180714_HARs.bed"

if [[ ! -f "$HAR_BED" ]]; then
  echo "[*] Downloading HARs (hg38) ..."
  $DL "$HAR_GZ" "$HAR_URL"
  echo "[*] Unzipping HARs ..."
  gunzip -f "$HAR_GZ"
else
  echo "[*] Using existing $HAR_BED"
fi

# Sanity check HARs format (BED3+ assumed)
if ! head -n1 "$HAR_BED" | awk 'BEGIN{FS="\t"} {exit !(NF>=3)}'; then
  echo "ERROR: HAR BED seems malformed (needs >=3 columns)" >&2; exit 1
fi

tr -d '\r' < GSE180714_HARs.bed \
| awk 'BEGIN{OFS="\t"} $0!~/^(#|track|browser)/ && NF>0 { gsub(/[ ]+/, "\t"); print }' \
| awk 'BEGIN{FS=OFS="\t"}
       NF>=3 && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $2<$3 { print }' \
> GSE180714_HARs.clean.bed

# Optional: sort (bedtools likes sorted input)
sort -k1,1 -k2,2n GSE180714_HARs.clean.bed -o GSE180714_HARs.bed

# -------------------- download & build gene BED --------------------
GENE_BED="genes.hg38.${ANNOT}.bed"
TSS_BED="genes.hg38.${ANNOT}.TSS.bed"

if [[ "$ANNOT" == "gencode" ]]; then
  GTF_GZ="gencode.v44.annotation.gtf.gz"
  GTF="gencode.v44.annotation.gtf"
  if [[ ! -f "$GTF" ]]; then
    echo "[*] Downloading GENCODE v44 (hg38) GTF ..."
    $DL "$GTF_GZ" "$GENCODE_URL"
    echo "[*] Unzipping GTF ..."
    gunzip -f "$GTF_GZ"
  else
    echo "[*] Using existing $GTF"
  fi

  if [[ ! -f "$GENE_BED" ]]; then
    echo "[*] Building gene BED from GENCODE GTF ..."
    awk 'BEGIN{FS=OFS="\t"}
         $3=="gene" {
           # Fields: chrom start end gene_id gene_name strand
           match($0,/gene_id "([^"]+)"/,a);
           match($0,/gene_name "([^"]+)"/,b);
           if(a[1]!=""){
             print $1,$4-1,$5,a[1],(b[1]!=""?b[1]:"."),$7
           }
         }' "$GTF" > "$GENE_BED"
  fi

  if [[ ! -f "$TSS_BED" ]]; then
    echo "[*] Computing TSS BED from gene BED ..."
    awk 'BEGIN{FS=OFS="\t"}
         ($6=="+"){print $1,$2,$2+1,$4,$5,$6}
         ($6=="-"){print $1,$3-1,$3,$4,$5,$6}' "$GENE_BED" > "$TSS_BED"
  fi

elif [[ "$ANNOT" == "refseq" ]]; then
  R_GZ="refGene.txt.gz"
  R_TXT="refGene.txt"
  if [[ ! -f "$R_TXT" ]]; then
    echo "[*] Downloading UCSC refGene table (hg38) ..."
    $DL "$R_GZ" "$REFSEQ_URL"
    echo "[*] Unzipping refGene ..."
    gunzip -f "$R_GZ"
  else
    echo "[*] Using existing $R_TXT"
  fi

  # refGene.txt fields (tab-separated):
  # 1:bin 2:name(refseq_id) 3:chrom 4:strand 5:txStart 6:txEnd ... 13:name2(geneSymbol)
  if [[ ! -f "$GENE_BED" ]]; then
    echo "[*] Building gene BED from refGene.txt ..."
    awk 'BEGIN{FS=OFS="\t"}
         {chrom=$3; strand=$4; start=$5; end=$6; id=$2; sym=$13;
          if(chrom ~ /^chr/){
            print chrom,start,end,id,(sym!=""?sym:"."),strand
          }}' "$R_TXT" \
      | sort -k1,1 -k2,2n \
      | awk 'BEGIN{FS=OFS="\t"}
             # Collapse to gene-level span per gene symbol if possible
             { key = $5"|"$6; if(!(key in minS) || $2<minS[key]) minS[key]=$2;
                               if(!(key in maxE) || $3>maxE[key]) maxE[key]=$3;
                               chr[key]=$1; id[key]=$4; sym[key]=$5; strand[key]=$6 }
             END{ for(k in chr) print chr[k],minS[k],maxE[k],id[k],sym[k],strand[k] }' \
      > "$GENE_BED"
  fi

  if [[ ! -f "$TSS_BED" ]]; then
    echo "[*] Computing TSS BED from refSeq gene BED ..."
    awk 'BEGIN{FS=OFS="\t"}
         ($6=="+"){print $1,$2,$2+1,$4,$5,$6}
         ($6=="-"){print $1,$3-1,$3,$4,$5,$6}' "$GENE_BED" > "$TSS_BED"
  fi
fi

# -------------------- intersections --------------------
if [[ "$MODE" == "overlap" || "$MODE" == "both" ]]; then
  OUT_BED="${OUTPREFIX}.genes_overlapping_HARs.bed"
  echo "[*] Intersecting genes with HARs (overlap) ..."
  bedtools intersect -u -a "$GENE_BED" -b "$HAR_BED" > "$OUT_BED"
  echo "[*] Wrote: $OUT_BED"

  if [[ "$ANNOT" == "gencode" ]]; then
    OUT_IDS="${OUTPREFIX}.ensembl_gene_ids.txt"
    awk 'BEGIN{FS=OFS="\t"} {print $4}' "$OUT_BED" | sort -u > "$OUT_IDS"
    echo "[*] Wrote: $OUT_IDS (unique Ensembl gene IDs)"
  else
    OUT_IDS="${OUTPREFIX}.refseq_ids_and_symbols.tsv"
    awk 'BEGIN{FS=OFS="\t"} {print $4,$5}' "$OUT_BED" | sort -u > "$OUT_IDS"
    echo "[*] Wrote: $OUT_IDS (RefSeq ID <tab> gene symbol)"
  fi
fi

if [[ "$MODE" == "nearest" || "$MODE" == "both" ]]; then
  OUT_NEAR="${OUTPREFIX}.nearest_TSS.tsv"
  echo "[*] Mapping each HAR to nearest TSS ..."
  # Output: HAR fields + nearest gene fields + distance
  bedtools closest -D a -a "$HAR_BED" -b "$TSS_BED" > "$OUT_NEAR"
  echo "[*] Wrote: $OUT_NEAR"
fi

echo "[done]"
