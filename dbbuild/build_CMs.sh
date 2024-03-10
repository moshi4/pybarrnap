#!/bin/bash

RFAM_VERSION="14.10"
DB_OUTDIR="./db/cm"

# Download RFAM CM
RFAMDB="Rfam.cm"
RFAMDB_URL="ftp://ftp.ebi.ac.uk/pub/databases/Rfam/${RFAM_VERSION}/${RFAMDB}.gz"
if [ ! -r "$RFAMDB" ]; then
    echo Downloading: $RFAMDB
    wget $RFAMDB_URL -q --show-progress
    gunzip ${RFAMDB}.gz
else
    echo Using existing file: $RFAMDB
fi

# Retrieve target accession rRNA CM
echo Retrieve target accession CM from CMs database...

echo Indexing $RFAMDB
cmfetch --index $RFAMDB

echo "Fetch Bacteria rRNA CMs (5S=RF00001, 16S=RF00177, 23S=RF02541)"
cmfetch $RFAMDB RF00001 >5S.bac.cm
cmfetch $RFAMDB RF00177 >16S.bac.cm
cmfetch $RFAMDB RF02541 >23S.bac.cm
BAC_CM_FILE="bac.cm"
cat *.bac.cm >$BAC_CM_FILE
sed -i -e "s/SSU_rRNA_bacteria/16S_rRNA/g" $BAC_CM_FILE
sed -i -e "s/LSU_rRNA_bacteria/23S_rRNA/g" $BAC_CM_FILE

echo "Fetch Archaea rRNA CMs (5S=RF00001, 5.8S=RF00002, 16S=RF01959, 23S=RF02540)"
cmfetch $RFAMDB RF00001 >5S.arc.cm
cmfetch $RFAMDB RF00002 >5_8S.arc.cm
cmfetch $RFAMDB RF01959 >16S.arc.cm
cmfetch $RFAMDB RF02540 >23S.arc.cm
ARC_CM_FILE="arc.cm"
cat *.arc.cm >$ARC_CM_FILE
sed -i -e "s/SSU_rRNA_archaea/16S_rRNA/g" $ARC_CM_FILE
sed -i -e "s/LSU_rRNA_archaea/23S_rRNA/g" $ARC_CM_FILE

echo "Fetch Eukaryote rRNA CMs (5S=RF00001, 5.8S=RF00002, 18S=RF01960, 28S=RF02543)"
cmfetch $RFAMDB RF00001 >5S.euk.cm
cmfetch $RFAMDB RF00002 >5_8S.euk.cm
cmfetch $RFAMDB RF01960 >18S.euk.cm
cmfetch $RFAMDB RF02543 >28S.euk.cm
EUK_CM_FILE="euk.cm"
cat *.euk.cm >$EUK_CM_FILE
sed -i -e "s/SSU_rRNA_eukarya/18S_rRNA/g" $EUK_CM_FILE
sed -i -e "s/LSU_rRNA_eukarya/28S_rRNA/g" $EUK_CM_FILE

# Build CM database
mkdir -p $DB_OUTDIR
for CM_FILE in $BAC_CM_FILE $ARC_CM_FILE $EUK_CM_FILE; do
    mv $CM_FILE $DB_OUTDIR
    cmpress -F ${DB_OUTDIR}/${CM_FILE}
done

# Remove unnecessary files
rm -f ${RFAMDB}.ssi
for KINGDOM in arc bac euk; do
    rm -f *${KINGDOM}.cm ${DB_OUTDIR}/${KINGDOM}.cm
done

# Show CM database files
echo -e "\nFinished building CM database:"
ls -1 ${DB_OUTDIR}/*.cm.*