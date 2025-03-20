#!/usr/bin/env sh

mkdir data/reference/Homo_sapiens
curl https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz -o data/reference/Homo_sapiens/chr21.fasta.gz
curl https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz -o data/reference/Homo_sapiens/chr22.fasta.gz
