# RNA Editing Datasets

This is a repository for downloading and simulating sequence datasets for studying RNA editing.

0. Create and activate Conda environment.

```sh
conda create env -f environment.yml
conda activate RNAEditingDatasets
```

1. Download sequences and annotations of *Homo sapiens* chromosomes 21 and 22.

```sh
sh src/download.sh
```

2. Simulate a diploid genome and transcriptome with SNPs and A-to-I transcript editions.

```sh
julia src/SNPpers.jl -z -o data/simulated data/reference/Homo_sapiens/chr21.fasta.gz
```

3. Simulate RNAseq reads.

```sh
# Coming soon...
```

---

![MiVEGEC](assets/MIVEGEC.png)
![IRD](assets/IRD.png)
