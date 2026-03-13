#!/usr/bin/env bash
# download_dbs.sh
#
# Downloads the three reference databases required by MetaFuncDecoder
# into the data/ directory.
#
# Usage:
#   bash download_dbs.sh
#   bash download_dbs.sh --data-dir /path/to/custom/dir
#
# Database sources:
#   Pfam2GO  — Gene Ontology Consortium (open access)
#              https://current.geneontology.org/ontology/external2go/pfam2go
#   CAZy     — dbCAN FamInfo (open access)
#              https://pro.unl.edu/dbCAN2/download_file.php?file=Databases/FamInfo.txt.04232020.xls
#   KEGG     — KOfam ko_list (KEGG license applies; free for academic use)
#              https://www.genome.jp/ftp/db/kofam/
#              NOTE: Download requires agreeing to KEGG's terms of use.
#              The ko_list file is downloaded here from the KOfam FTP.
#              For commercial use, consult https://www.kegg.jp/kegg/legal.html

set -euo pipefail

# --- defaults -------------------------------------------------------------

DATA_DIR="data"

# --- argument parsing -----------------------------------------------------

while [[ $# -gt 0 ]]; do
    case "$1" in
        --data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        -h|--help)
            grep "^#" "$0" | sed 's/^# \{0,2\}//'
            exit 0
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

# --- setup ----------------------------------------------------------------

mkdir -p "$DATA_DIR"
echo "Downloading reference databases to: $DATA_DIR/"
echo ""

# --- Pfam2GO --------------------------------------------------------------

PFAM2GO="$DATA_DIR/pfam2go.txt"
echo "[1/3] Pfam2GO mapping"
echo "  Source: https://current.geneontology.org/ontology/external2go/pfam2go"
if [[ -f "$PFAM2GO" ]]; then
    echo "  Skipping — $PFAM2GO already exists (delete to re-download)"
else
    curl -fsSL \
        "https://current.geneontology.org/ontology/external2go/pfam2go" \
        -o "$PFAM2GO"
    echo "  Saved: $PFAM2GO"
fi
echo ""

# --- CAZy FamInfo ---------------------------------------------------------

CAZY_XLS="$DATA_DIR/FamInfo.txt.04232020.xls"
CAZY_INFO="$DATA_DIR/FamInfo.txt.04232020.tsv"
echo "[2/3] CAZy family information (FamInfo)"
echo "  Source: https://pro.unl.edu/dbCAN2/download_file.php?file=Databases/FamInfo.txt.04232020.xls"
if [[ -f "$CAZY_INFO" ]]; then
    echo "  Skipping — $CAZY_INFO already exists (delete to re-download)"
else
    curl -fsSL \
        "https://pro.unl.edu/dbCAN2/download_file.php?file=Databases/FamInfo.txt.04232020.xls" \
        -o "$CAZY_XLS"
    python -c "
import pandas as pd, sys
try:
    df = pd.read_excel('$CAZY_XLS', sheet_name=0)
    df.to_csv('$CAZY_INFO', sep='\t', index=False)
    print('  Converted to TSV: $CAZY_INFO')
except Exception as e:
    print('  Error converting XLS to TSV:', e)
    print('  Try: pip install xlrd')
    sys.exit(1)
"
    rm -f "$CAZY_XLS"
fi
echo ""

# --- KEGG KOfam ko_list ---------------------------------------------------

KO_LIST="$DATA_DIR/ko_list"
KO_GZ="$DATA_DIR/ko_list.gz"
echo "[3/3] KEGG KOfam ko_list"
echo "  Source: https://www.genome.jp/ftp/db/kofam/ko_list.gz"
echo "  NOTE: KEGG data is free for academic use."
echo "  For commercial use see: https://www.kegg.jp/kegg/legal.html"
if [[ -f "$KO_LIST" ]]; then
    echo "  Skipping — $KO_LIST already exists (delete to re-download)"
else
    curl -fsSL \
        "https://www.genome.jp/ftp/db/kofam/ko_list.gz" \
        -o "$KO_GZ"
    gunzip "$KO_GZ"
    echo "  Saved: $KO_LIST"
fi
echo ""

# --- summary --------------------------------------------------------------

echo "Done. Reference databases:"
for f in "$PFAM2GO" "$CAZY_INFO" "$KO_LIST"; do
    if [[ -f "$f" ]]; then
        size=$(du -sh "$f" | cut -f1)
        echo "  $size  $f"
    else
        echo "  MISSING  $f"
    fi
done
echo ""
echo "Run the tool:"
echo "  python metafuncdecoder.py --ko-file your.assembled.ko --function carbon_cycling"
