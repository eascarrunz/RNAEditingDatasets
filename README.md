# RNA Editing Datasets

This is a repository for downloading and simulating sequence datasets for studying RNA editing.

## ssnpper

**ssnpper** (/snɪ́pə/) performs **S**imulation of **SNP**s **P**lus **E**ditions of **R**NA.

This script takes in sequences from a reference genome and generates a new genome with point mutations (SNP) and a transcriptome with RNA editions.
The corresponding output files are tagged "genome" and "transcriptome".

By default, the script simulates a diploid genome and transcriptome from an unphased reference genome. Output files for each simulated haplotype are tagged "hap" with an ID number.
The ploidy can be adjusted with the `--ploidy` option; the ploidy should be set to 1 for simulations on phased reference genomes.

SNPs are simulated only with the JC96 model, and only A-to-I RNA editions are implemented.

### Setup

#### With local installation of Julia

The script was developed with Julia 1.10.9. It is recommended that you install juliaup to be able to use choose the precise version of Julia that you wish to use on your machine. This is done by default with the officially recommended Julia installation method:

```sh
curl -fsSL https://install.julialang.org | sh
```

Then you can download and switch to version 1.10.9

```sh
juliaup add 1.10.9
juliaup default 1.10.9
```

After those steps, it suffices to clone the repository. The script will automatically download and install the dependencies to a dedicated environment based on the Project.toml file.

#### With Docker

A Dockerfile is provided. To build it, use the following command:

```sh
docker build -t ssnpper -f docker/ssnpper/Dockerfile .
```

The image can be run like so:

```sh
docker run -v .:/work -it ssnpper
```

Within the Docker container, use the `ssnpper` alias to run the script instead of a Julia command. For example:

```sh
ssnpper -z -o data/simulated data/reference/Homo_sapiens/chr21.fasta.gz
```

### Usage

1. Download sequences and annotations of *Homo sapiens* chromosomes 21 and 22.

```sh
sh src/download.sh
```

2. Simulate a diploid genome and transcriptome with SNPs and A-to-I transcript editions.

```sh
julia ssnpper.jl -z -o data/simulated data/reference/Homo_sapiens/chr21.fasta.gz
```

3. Simulate RNAseq reads.

```sh
# Coming soon...
```

---

![MiVEGEC](assets/MIVEGEC.png)
![IRD](assets/IRD.png)
