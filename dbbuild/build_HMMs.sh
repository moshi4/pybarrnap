#!/bin/bash

RFAM_VERSION="14.10"
DB_OUTDIR="./db/hmm"

# Download RFAM annotated seed alignments
RFAMDB="Rfam.seed"
RFAMDB_URL="ftp://ftp.ebi.ac.uk/pub/databases/Rfam/${RFAM_VERSION}/${RFAMDB}.gz"
if [ ! -r "$RFAMDB" ]; then
    echo Downloading: $RFAMDB
    wget $RFAMDB_URL -q --show-progress
    gunzip ${RFAMDB}.gz
else
    echo Using existing file: $RFAMDB
fi

# Retrieve target accession rRNA MSA
echo Retrieve target accession MSA from MSAs database...

echo Indexing $RFAMDB
esl-afetch --index $RFAMDB

echo "Fetch Bacteria rRNA MSAs (5S=RF00001, 16S=RF00177, 23S=RF02541)"
esl-afetch $RFAMDB RF00001 >5S.bac.aln
esl-afetch $RFAMDB RF00177 >16S.bac.aln
esl-afetch $RFAMDB RF02541 >23S.bac.aln

echo "Fetch Archaea rRNA MSAs (5S=RF00001, 5.8S=RF00002, 16S=RF01959, 23S=RF02540)"
esl-afetch $RFAMDB RF00001 >5S.arc.aln
esl-afetch $RFAMDB RF00002 >5_8S.arc.aln
esl-afetch $RFAMDB RF01959 >16S.arc.aln
esl-afetch $RFAMDB RF02540 >23S.arc.aln

echo "Fetch Eukaryote rRNA MSAs (5S=RF00001, 5.8S=RF00002, 18S=RF01960, 28S=RF02543)"
esl-afetch $RFAMDB RF00001 >5S.euk.aln
esl-afetch $RFAMDB RF00002 >5_8S.euk.aln
esl-afetch $RFAMDB RF01960 >18S.euk.aln
esl-afetch $RFAMDB RF02543 >28S.euk.aln

# Build HMM database
mkdir -p $DB_OUTDIR
for KINGDOM in bac arc euk; do
    for TYPE in 5S 5_8S 16S 23S 18S 28S; do
        ID="${TYPE}.${KINGDOM}"
        if [ -s "${ID}.aln" ]; then
            echo "*** $ID ***"
            hmmbuild --rna -n ${TYPE}_rRNA ${ID}.hmm ${ID}.aln
        fi
    done
    cat *.${KINGDOM}.hmm >${DB_OUTDIR}/${KINGDOM}.hmm
done

# Remove unnecessary files
rm -f ${RFAMDB}.ssi
for KINGDOM in arc bac euk; do
    rm -f *.${KINGDOM}.aln *.${KINGDOM}.hmm
done

# Show HMM database files
echo -e "\nFinished building HMM profiles:"
ls -1 ${DB_OUTDIR}/*.hmm
