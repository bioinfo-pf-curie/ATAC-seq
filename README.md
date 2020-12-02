# ATAC-seq

**Institut Curie - Nextflow ATAC-seq analysis pipeline**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)
[![MultiQC](https://img.shields.io/badge/MultiQC-1.8-blue.svg)](https://multiqc.info/)
[![Install with](https://anaconda.org/anaconda/conda-build/badges/installer/conda.svg)](https://conda.anaconda.org/anaconda)
[![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)](https://singularity.lbl.gov/)
[![Docker Container available](https://img.shields.io/badge/docker-available-003399.svg)](https://www.docker.com/)


**UNDER CONSTRUCTION - DO NOT USE**

### Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. 
It comes with containers making installation trivial and results highly reproducible.
The current workflow was initiated from the [nf-core ATAC-seq pipeline](https://github.com/nf-core/chipseq). See the nf-core project from details on [guidelines](https://nf-co.re/).

### Pipeline Summary

1. Run quality control of raw sequencing reads ([`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Align reads on reference genome ([`BWA`](http://bio-bwa.sourceforge.net/) / [`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) / [`STAR`](https://github.com/alexdobin/STAR))
3. Sort aligned reads ([`SAMTools`](http://www.htslib.org/))
4. Mark duplicates ([`Picard`](https://broadinstitute.github.io/picard/))
5. Library complexity analysis ([`Preseq`](http://smithlabresearch.org/software/preseq/))
6. Filtering aligned BAM files ([`SAMTools`](http://www.htslib.org/) & [`BAMTools`](https://github.com/pezmaster31/bamtools))
   - reads mapped to mitochondrial DNA
   - reads mapped to blacklisted regions
   - reads marked as duplicates
   - reads that arent marked as primary alignments
   - reads that are unmapped
   - reads mapped with a low mapping quality (multiple hits, secondary alignments, etc.)
7. Create normalized bigWig file ([`deepTools`](https://deeptools.readthedocs.io/en/develop/index.html))
8. Peak calling ([`MACS2`](https://github.com/taoliu/MACS))
9. Peak annotation and QC ([`HOMER`](http://homer.ucsd.edu/homer/ngs/annotation.html))
10. Results summary ([`MultiQC`](https://multiqc.info/))

### Quick help

```bash

N E X T F L O W  ~  version 20.01.0
======================================================================
ATAC-seq v1.0.1
======================================================================

Usage:

nextflow run main.nf --reads '*_R{1,2}.fastq.gz' -profile conda --genomeAnnotationPath '/data/annotations/pipelines' --genome 'hg19'
nextflow run main.nf --samplePlan 'sample_plan.csv' --design 'design.csv' -profile conda --genomeAnnotationPath '/data/annotations/pipelines' --genome 'hg19'

Mandatory arguments:
--reads [file]                     Path to input data (must be surrounded with quotes)
--samplePlan [file]                Path to sample plan file if '--reads' is not specified
--genome [str]                     Name of genome reference. See the `--genomeAnnotationPath` to defined the annotations path
-profile [str]                     Configuration profile to use. Can use multiple (comma separated)

Inputs:
--singleEnd [bool]                 Specifies that the input is single end reads
--fragmentSize [int]               Estimated fragment length used to extend single-end reads. Default: 200

References           If not specified in the configuration file or you wish to overwrite any of the references given by the --genome field
--fasta [file]                     Path to Fasta reference

Alignment:
--aligner [str]                    Alignment tool to use ['bwa-mem', 'star', 'bowtie2']. Default: 'bwa-mem'
--saveAlignedIntermediates [bool]  Save all intermediates mapping files. Default: false  
--starIndex [file]                 Index for STAR aligner
--bwaIndex [file]                  Index for Bwa-mem aligner
--bowtie2Index [file]              Index for Bowtie2 aligner

Filtering:
--mapq [int]                       Minimum mapping quality to consider. Default: false
--keepDups [bool]                  Do not remove duplicates afer marking. Default: false
--keepSingleton [bool]             Keep unpaired reads. Default: false
--blacklist [file]                 Path to black list regions (.bed).

Annotation:          If not specified in the configuration file or you wish to overwrite any of the references given by the --genome field
--geneBed [file]                   BED annotation file with gene coordinate.
--gtf [file]                       GTF annotation file. Used in HOMER peak annotation
--effGenomeSize [int]              Effective Genome size
--tssSize [int]                    Distance (upstream/downstream) to transcription start point to consider. Default: 2000

Skip options:        All are false by default
--skipFastqc [bool]                Skips fastQC
--skipPreseq [bool]                Skips preseq QC
--skipDeepTools [bool]             Skips deeptools QC
--skipPeakcalling [bool]           Skips peak calling
--skipPeakanno [bool]              Skips peak annotation
--skipFeatCounts [bool]            Skips feature count
--skipMultiQC [bool]               Skips MultiQC step

Other options:
--outdir [file]                    The output directory where the results will be saved
--email [str]                      Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
-name [str]                        Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

=======================================================
Available Profiles
  -profile test                    Run the test dataset
  -profile conda                   Build a new conda environment before running the pipeline. Use `--condaCacheDir` to define the conda cache path
  -profile multiconda              Build a new conda environment per process before running the pipeline. Use `--condaCacheDir` to define the conda cache path
  -profile path                    Use the installation path defined for all tools. Use `--globalPath` to define the insallation path
  -profile multipath               Use the installation paths defined for each tool. Use `--globalPath` to define the insallation path
  -profile docker                  Use the Docker images for each process
  -profile singularity             Use the Singularity images for each process. Use `--singularityPath` to define the insallation path
  -profile cluster                 Run the workflow on the cluster, instead of locally

```

### Quick run

The pipeline can be run on any infrastructure from a list of input files or from a sample plan as follow

#### Run the pipeline on a test dataset
See the conf/test.conf to set your test dataset.

```
nextflow run main.nf -profile test,conda

```

#### Run the pipeline from a `sample plan` and a `design` file
```
nextflow run main.nf --samplePlan MY_SAMPLE_PLAN --design MY_DESIGN --genome 'hg19' --genomeAnnotationPath ANNOTATION_PATH --outdir MY_OUTPUT_DIR

```

### Defining the '-profile'

By default (whithout any profile), Nextflow will excute the pipeline locally, expecting that all tools are available from your `PATH` variable.

In addition, we set up a few profiles that should allow you i/ to use containers instead of local installation, ii/ to run the pipeline on a cluster instead of on a local architecture.
The description of each profile is available on the help message (see above).

Here are a few examples of how to set the profile option.

```
## Run the pipeline locally, using a global environment where all tools are installed (build by conda for instance)
-profile path --globalPath INSTALLATION_PATH

## Run the pipeline on the cluster, using the Singularity containers
-profile cluster,singularity --singularityPath SINGULARITY_PATH

## Run the pipeline on the cluster, building a new conda environment
-profile cluster,conda --condaCacheDir CONDA_CACHE

```

### Sample Plan

A sample plan is a csv file (comma separated) that list all samples with their biological IDs.
The sample plan is expected to be created as below :

SAMPLE_ID | SAMPLE_NAME | FASTQ_R1 [Path to R1.fastq file] | FASTQ_R2 [For paired end, path to Read 2 fastq]

### Full Documentation

1. [Installation](docs/installation.md)
2. [Reference genomes](docs/reference_genomes.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

#### Credits

This pipeline has been written by the bioinformatics platform of the Institut Curie (Valentin Laroche, Nicolas Servant)

#### Contacts

For any question, bug or suggestion, please use the issues system or contact the bioinformatics core facility.
