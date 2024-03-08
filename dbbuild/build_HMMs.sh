#!/bin/bash

CPUS=$(grep -c bogomips /proc/cpuinfo)
DB_OUTDIR="./db"

# Download RFAM annotated seed alignments
RFAM="Rfam.seed"
RFAMURL="ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/${RFAM}.gz"
if [ ! -r "$RFAM" ]; then
  echo Downloading: $RFAM
  wget $RFAMURL -q --show-progress
  gunzip ${RFAM}.gz
else
  echo Using existing file: $RFAM
fi

# Retrieve target accession rRNA MSA
echo Retrieve target accession MSA from MSA database...

echo Indexing $RFAM
esl-afetch --index $RFAM

echo "Bacteria (5S=RF00001, 16S=RF00177, 23S=RF02541)"
esl-afetch $RFAM RF00001 >5S.bac.aln
esl-afetch $RFAM RF00177 >16S.bac.aln
esl-afetch $RFAM RF02541 >23S.bac.aln

echo "Archaea (5S=RF00001, 5.8S=RF00002, 16S=RF01959, 23S=RF02540)"
esl-afetch $RFAM RF00001 >5S.arc.aln
esl-afetch $RFAM RF00002 >5_8S.arc.aln
esl-afetch $RFAM RF01959 >16S.arc.aln
esl-afetch $RFAM RF02540 >23S.arc.aln

echo "Eukaryote (5S=RF00001, 5.8S=RF00002, 18S=RF01960, 28S=RF02543)"
esl-afetch $RFAM RF00001 >5S.euk.aln
esl-afetch $RFAM RF00002 >5_8S.euk.aln
esl-afetch $RFAM RF01960 >18S.euk.aln
esl-afetch $RFAM RF02543 >28S.euk.aln

# Check metazoan mitochondria alignment file exists
FILE="12S.mito.aln"
if [ ! -r "$FILE" ]; then
  echo "Missing included $FILE file."
  exit 1
fi
FILE="16S.mito.aln"
if [ ! -r "$FILE" ]; then
  echo "Missing included $FILE file."
  exit 1
fi

# Build HMM profile
mkdir -p $DB_OUTDIR
for KINGDOM in arc bac euk mito; do
  for TYPE in 5S 5_8S 12S 16S 23S 18S 28S; do
    ID="${TYPE}.${KINGDOM}"
    if [ -s "${ID}.aln" ]; then
      echo "*** $ID ***"
      hmmbuild --cpu $CPUS --rna -n ${TYPE}_rRNA ${ID}.hmm ${ID}.aln
    fi
  done
  cat *.${KINGDOM}.hmm >${DB_OUTDIR}/${KINGDOM}.hmm
done

# Remove unnecessary files
rm -f ${RFAM}.ssi
for KINGDOM in arc bac euk; do
  rm -f *.${KINGDOM}.aln *.${KINGDOM}.hmm
done
rm -f *.mito.hmm

# Show HMM profile files
echo -e "\nFinished building HMM profiles:"
ls -1 ${DB_OUTDIR}/*.hmm
