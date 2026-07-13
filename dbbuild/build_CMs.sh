#!/bin/bash

RFAM_VERSION="15.1"
DB_OUTDIR="./db/cm"
mkdir -p $DB_OUTDIR

# Download RFAM CM and clan file
RFAMDB="Rfam.cm"
RFAM_CLANIN="Rfam.clanin"
RFAMDB_URL="ftp://ftp.ebi.ac.uk/pub/databases/Rfam/${RFAM_VERSION}/${RFAMDB}.gz"
RFAM_CLANIN_URL="ftp://ftp.ebi.ac.uk/pub/databases/Rfam/${RFAM_VERSION}/${RFAM_CLANIN}"
if [ ! -r "$RFAMDB" ]; then
    echo Downloading: $RFAMDB
    wget $RFAMDB_URL -q --show-progress
    gunzip ${RFAMDB}.gz
else
    echo Using existing file: $RFAMDB
fi
if [ ! -r "$RFAM_CLANIN" ]; then
    echo Downloading: $RFAM_CLANIN
    wget $RFAM_CLANIN_URL -q --show-progress
else
    echo Using existing file: $RFAM_CLANIN
fi

# Define Accessions

    # Previosly, cmfetch was used with accessions. Since the clanin file works with nams, but the names are the same as the ones that
    # can be used by cmfetch, we define the names here once and use them for both cm and clanin.

# (5S=RF00001=5S_rRNA, 16S=RF00177=SSU_rRNA_bacteria, 23S=RF02541=LSU_rRNA_bacteria)
BAC_MODELS="5S_rRNA SSU_rRNA_bacteria LSU_rRNA_bacteria"
# (5S=RF00001=5S_rRNA, 5.8S=RF00002=5_8S_rRNA, 16S=RF01959=SSU_rRNA_archaea, 23S=RF02540=LSU_rRNA_archaea)
ARC_MODELS="5S_rRNA 5_8S_rRNA SSU_rRNA_archaea LSU_rRNA_archaea"
# (5S=RF00001=5S_rRNA, 5.8S=RF00002=5_8S_rRNA, 18S=RF01960=SSU_rRNA_eukarya, 28S=RF02543=LSU_rRNA_eukarya)
EUK_MODELS="5S_rRNA 5_8S_rRNA SSU_rRNA_eukarya LSU_rRNA_eukarya"
# Deduplicate shared models for the all-kingdom database.
ALL_MODELS="$(
    printf '%s\n' $BAC_MODELS $ARC_MODELS $EUK_MODELS |
        sort -u |
        tr '\n' ' '
)"


# Retrieve target accession rRNA CM
echo Retrieve target accession CM from CMs database...

echo Indexing $RFAMDB
cmfetch --index $RFAMDB

fetch_models() {
    local models="$1"
    local outfile="$2"

    for model in $models; do
        cmfetch "$RFAMDB" "$model"
    done >"$outfile"
}

echo "Fetch Bacteria rRNA CMs (5S=RF00001, 16S=RF00177, 23S=RF02541)"
BAC_CM_FILE="bac.cm"
fetch_models "$BAC_MODELS" "$BAC_CM_FILE"

echo "Fetch Archaea rRNA CMs (5S=RF00001, 5.8S=RF00002, 16S=RF01959, 23S=RF02540)"
ARC_CM_FILE="arc.cm"
fetch_models "$ARC_MODELS" "$ARC_CM_FILE"

echo "Fetch Eukaryote rRNA CMs (5S=RF00001, 5.8S=RF00002, 18S=RF01960, 28S=RF02543)"
EUK_CM_FILE="euk.cm"
fetch_models "$EUK_MODELS" "$EUK_CM_FILE"

echo "Fetch All Kingdoms rRNA CMs (RF00001, RF00002, RF00177, RF01959, RF01960, RF02540, RF02541, RF02543)"
ALL_CM_FILE="all.cm"
fetch_models "$ALL_MODELS" "$ALL_CM_FILE"

# Filter clanin file

    # The clanin file looks like this:
    # CL00001	alpha_tmRNA	beta_tmRNA	cyano_tmRNA	mt-tmRNA	tmRNA	tRNA	tRNA-Sec
    # CL00002	RNaseP-T	RNaseP_arch	RNaseP_bact_a	RNaseP_bact_b	RNaseP_nuc	RNase_MRP	RNase_P

    # So a space-separated file where each row is one clan, the first column is the ID, the other columns are names of identical RNAs.

    # BEGIN -> run once before looking at any file lines.
    # split(models, model_arr, " ") turn the space-separated list into an array.
    # for loop builds a lookup table, marked as "1" -> keep this

    # for every column beginning with 2
    # if the column is in the lookup table, add the whole line to output
    # if count > 0, print the line

filter_clanin() {
    local models="$1"
    local outfile="$2"
    awk -v models="$models" '
    BEGIN {
        split(models, model_arr, " ")
        for (idx in model_arr) {
            keep[model_arr[idx]] = 1
        }
    }
    {
        line = $1
        count = 0
        for (idx = 2; idx <= NF; idx++) {
            if ($idx in keep) {
                line = line "\t" $idx
                count++
            }
        }
        if (count > 0) {
            print line
        }
    }
' "$RFAM_CLANIN" >"$outfile"
}
echo "Fetch Bacteria clanins"
filter_clanin "$BAC_MODELS" "${DB_OUTDIR}/bac.clanin"
echo "Fetch Archaeal clanins"
filter_clanin "$ARC_MODELS" "${DB_OUTDIR}/arc.clanin"
echo "Fetch Eukaryote clanins"
filter_clanin "$EUK_MODELS" "${DB_OUTDIR}/euk.clanin"
echo "Fetch All Kingdoms clanins"
filter_clanin "$ALL_MODELS" "${DB_OUTDIR}/all.clanin"

# Build CM database
for CM_FILE in $BAC_CM_FILE $ARC_CM_FILE $EUK_CM_FILE $ALL_CM_FILE; do
    mv $CM_FILE $DB_OUTDIR
    cmpress -F ${DB_OUTDIR}/${CM_FILE}
done

# Remove unnecessary files
rm -f ${RFAMDB}.ssi
for KINGDOM in arc bac euk all; do
    rm -f *${KINGDOM}.cm ${DB_OUTDIR}/${KINGDOM}.cm
done

# Show CM database files
echo -e "\nFinished building CM database:"
ls -1 ${DB_OUTDIR}/*.cm.* ${DB_OUTDIR}/*.clanin
