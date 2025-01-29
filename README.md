# SDRranger
A tool for processing SDR-seq (targeted Single-cell DNA and RNA sequencing) fastq files into annotated BAM files and count matrices.

[![DOI](https://zenodo.org/badge/716143278.svg)](https://doi.org/10.5281/zenodo.14762617)

<p align="center">
  <img src="doc/SDRLoneRanger.png" width=50% alt="The Single-Cell DNA RNA Lone Ranger">
</p>
<p align="center">The Single-Cell DNA RNA Lone Ranger</p>

## Installation

SDRranger works in Linux, and has been tested on el8 and CentOS 7. 

You can install SDRranger from github directly into your local (conda or virtual) environment using pip:
```
pip install git+https://github.com/hawkjo/SDRranger.git
```
This typically takes a few minutes. 

[STAR aligner](https://github.com/alexdobin/STAR)  needs to be installed independently as well.

## Usage

The basic usage for SDRranger can be displayed at any time via `SDRranger --help`:
```
Usage:
  SDRranger count_RNA        <fastq_dir> (--STAR-ref-dir=<> | --STAR-output=<>...) --config=<> [--output-dir=<>] [--threads=<>] [-v | -vv | -vvv]
  SDRranger count_gDNA       <fastq_dir> (--STAR-ref-dir=<> | --STAR-output=<>...) --config=<> [--output-dir=<>] [--threads=<>] [-v | -vv | -vvv]
  SDRranger preprocess_gDNA  <fastq_dir> --config=<> [--output-dir=<>] [--threads=<>] [-v | -vv | -vvv]
  SDRranger count_matrix       <SDR_bam_file> --output-dir=<> [--threads=<>] [-v | -vv | -vvv]
  SDRranger simulate_reads   --config=<> --fastq-prefix=<> --nreads=<> [--unique-umis=<>] [--seed=<>] [--error-probability=<>] [--substitution-probability=<>] [--insertion-probability=<>] [-v | -vv | -vvv]

Options:
  --STAR-ref-dir=<>:              Path to directory with STAR index.
  --STAR-output=<>:               Path to STAR output file (BAM/SAM). Can be repeated multiple times,
                                    in which case the order must correspond to the lexicographic ordering
                                    of paired FASTQ files in <fastq_dir>.
  --config=<>:                    Path to JSON configuration.
  --output-dir=<>:                Path to output directory [default: .].
  --threads=<>:                   Number of threads [default: 1].
  -v:                             Verbose output.
  --fastq-prefix=<>:              Prefix for output FASTQ files.
  --nreads=<>:                    Number of reads to simulate.
  --unique-umis=<>:               Fraction of all reads that have unique UMIs [default: 0.5].
  --seed=<>:                      Random seed [default: 42].
  --error-probability=<>:         Probability of an error occurring [default: 0.1]. Set to a negative number to
                                    always introduce as many errors as allowed by the configuration.
  --substitution-probability=<>:  Probability of generating a substitution as opposed to an indel [default: 0.7].
  --insertion-probability=<>:     Probability of generating an insertion as opposed to a deletion when generating an indel [default: 0.5].
  -h --help                       Show this screen.
  --version                       Show version.

Commands:
  preprocess_gDNA  Preprocess Genomic gDNA files such that STAR can be run on the output
  count_gDNA       Process and count Genomic gDNA files
  count_RNA        Process and count Transcriptomic RNA files
  count_matrix     Build a count matrix (or matrices) from an existing bam file
  simulate_reads   Generate synthetic sequencing reads given a barcode configuration
```
As this shows, the typical workflow is performed on data separated into gDNA and RNA fastq files, as these two filetypes have different barcode structures and semantics, and require different handling.

The barcode details are input via a json configuration file, of which standard gDNA and RNA versions can be found in the `examples` folder. 

STAR references need to be prebuilt and their top directory input as a parameter.

### Outputs
The primary outputs from SDR ranger are:
* An annotated BAM file
* A read count matrix
* A UMI count matrix (for RNA)

#### BAM file tags
The BAM file is annotated with custom tags that have been created in the style of current community standards. These are:

| Tag | Meaning |
|--|--|
| CB | Cell barcode |
| CR | Raw, uncorrected cell barcode |
| UB | UMI |
| UR | Raw, uncorrected UMI |
| FL | Combined length of the linker sequences |

The cell barcode tag contains all pieces of the cell barcode, including the sample barcode, concatenated with periods.

## Examples

The `examples` folder contains small example gDNA and RNA datasets and corresponding example scripts and configuration files. The `cmd_gDNA_json.sh` and `cmd_cDNA_json.sh` files demonstrate proper syntax for their respective datasets and are runnable directly from within the examples folder. They each take about 20 seconds to run.
