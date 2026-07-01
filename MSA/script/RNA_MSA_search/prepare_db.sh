#!/bin/bash
# Build the local RNA MSA databases (Rfam + RNAcentral) used by msa_gen.py.
# Run this script from its own directory (MSA/script/RNA_MSA_search); it writes
# the databases to ./db there, which is exactly where msa_generation.py looks
# for them.
#
# Requirements on PATH: mmseqs, makehmmerdb, esl-sfetch (HMMER / easel), wget.
set -e

# --- Settings ---
DB_DIR="./db"
THREADS=$(nproc)
TMP_DIR="./tmp_mmseqs"

# Memory limit for mmseqs (use ~80% of the available RAM).
# e.g. 50G for a 64GB session, 100G for a 128GB session.
MEM_LIMIT="50G"

mkdir -p $DB_DIR
mkdir -p $TMP_DIR

echo "=== [1/2] Processing Rfam (latest) ==="
mkdir -p $DB_DIR/rfam
cd $DB_DIR/rfam

if [ ! -f "Rfam.fa.gz" ]; then
  wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/Rfam.fa.gz
fi
if [ ! -f "Rfam.fa" ]; then
  gunzip -k Rfam.fa.gz
fi

# Rfam is small, so no split option is needed.
if [ ! -f "rfam_clust_rep_seq.fasta" ]; then
  echo "Clustering Rfam..."
  mmseqs easy-cluster Rfam.fa rfam_clust $TMP_DIR \
    --min-seq-id 0.9 -c 0.8 --cov-mode 1 --threads $THREADS
fi

if [ ! -f "rfam_v_latest.mdf" ]; then
  echo "Building HMMER DB for Rfam..."
  makehmmerdb rfam_clust_rep_seq.fasta rfam_v_latest.mdf
fi

if [ ! -f "rfam_clust_rep_seq.fasta.ssi" ]; then
  esl-sfetch --index rfam_clust_rep_seq.fasta
fi
cd ../..

echo "=== [2/2] Processing RNAcentral (latest active) ==="
mkdir -p $DB_DIR/rnacentral
cd $DB_DIR/rnacentral

if [ ! -f "rnacentral_active.fasta.gz" ]; then
  wget ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_active.fasta.gz
fi
if [ ! -f "rnacentral_active.fasta" ]; then
  gunzip -k rnacentral_active.fasta.gz
fi

# RNAcentral is very large, so add --split-memory-limit.
if [ ! -f "rnacentral_clust_rep_seq.fasta" ]; then
  echo "Clustering RNAcentral (large dataset; split logic applied)..."
  mmseqs easy-linclust rnacentral_active.fasta rnacentral_clust $TMP_DIR \
    --min-seq-id 0.9 -c 0.8 --cov-mode 1 --threads $THREADS \
    --split-memory-limit $MEM_LIMIT
fi

if [ ! -f "rnacentral_v_latest.mdf" ]; then
  echo "Building HMMER DB for RNAcentral..."
  makehmmerdb rnacentral_clust_rep_seq.fasta rnacentral_v_latest.mdf
fi

if [ ! -f "rnacentral_clust_rep_seq.fasta.ssi" ]; then
  esl-sfetch --index rnacentral_clust_rep_seq.fasta
fi
cd ../..

echo "=== All databases prepared successfully ==="
