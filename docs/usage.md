# Usage

## Table of contents

* [Introduction](#general-nextflow-info)
* [Running the pipeline](#running-the-pipeline)
* [Main arguments](#main-arguments)
    * [`--reads`](#--reads)
    * [`--samplePlan`](#--samplePlan)
    * [`--singleEnd`](#--singleend)
    * [`--fragmentSize`](#--fragmentSize)
* [Reference genomes](#reference-genomes)
    * [`--genome`](#--genome)
    * [`--genomeAnnotationPath`](#--genomeAnnotationPath)
* [Alignment](#alignment)
    * [`--aligner`](#--aligner)
* [Filtering](#filtering)
    * [`--mapq`](#--mapq)
    * [`--keepDups`](#--keepDups)
	* [`--keepMito`](#--keepMito)
	* [`--keepSingleton`](#--keepSingleton)
* [Peak calling](#peak-calling)
    * [`--caller`](#--caller)
	* [`--tn5sites`](#--tn5sites)
	* [`--extsize`](#--extisze)
* [Profiles](#profiles)
* [Job resources](#job-resources)
* [Automatic resubmission](#automatic-resubmission)
* [Custom resource requests](#custom-resource-requests)
* [Other command line parameters](#other-command-line-parameters)
    * [`--skip*`](#--skip*)
    * [`--metadata`](#--metadta)
    * [`--outDir`](#--outDir)
    * [`-name`](#-name-single-dash)
    * [`-resume`](#-resume-single-dash)
    * [`-c`](#-c-single-dash)
    * [`--maxMemory`](#--maxMemory)
    * [`--maxTime`](#--maxTime)
    * [`--maxCpus`](#--maxCpus)
    * [`--multiqcConfig`](#--multiqcConfig)


## General Nextflow info
Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline
The typical command for running the pipeline is as follows:
```bash
nextflow run main.nf --reads '*_R{1,2}.fastq.gz' -profile 'singularity'
```

This will launch the pipeline with the `singularity` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

You can change the output director using the `--outDir/-w` options.

## Main arguments

### `--reads`
Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`


### `--samplePlan`
Use this to specify a sample plan file instead of a regular expression to find fastq files. For example :

```bash
--samplePlan 'path/to/data/sample_plan.csv
```

The sample plan is a csv file with the following information :

Sample ID | Sample Name | Path to R1 fastq file | Path to R2 fastq file

### `--singleEnd`

By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--singleEnd` on the command line when you launch the pipeline. A normal glob pattern, enclosed 
in quotation marks, can then be used for `--reads`. For example:

```bash
--singleEnd --reads '*.fastq.gz'
```

### `--fragmentSize`

For `--singleEnd` data only. Specify the expected fragment size in the experiments (default: 200).  
This information is used to extend the reads to fragment length in some analysis.  

```bash
--singleEnd --fragmentSize 200 --reads '*.fastq.gz'
```

## Reference Genomes

All information about genomes and annotation are available in `conf/genomes.config`.

### `-genome`

There are different species supported in the genomes references file. To run the pipeline, you must specify which to use with the `--genome` flag.

You can find the keys to specify the genomes in the [genomes config file](../conf/genomes.config). Common genomes that are supported are:

* Human
  * `--genome hg38`
* Mouse
  * `--genome mm10`

> There are numerous others - check the config file for more.

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the genomes resource. 
The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'hg19' {
      fasta         = '<path to the genome fasta file>'
      chrsize       = '<path to the chromosome lenght file>'
      bwaIndex      = '<path to the BWA-mem index files>'
      bowtie2Index  = '<path to the Bowtie2 index files>'
      geneBed       = '<path to a gene annotation file in bed format>'
      gtf           = '<path to annotation file in gtf format>'
      effGenomeSize = '<effective genome size>'
      blacklist     = '<path to black list region file>'
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

Note that these paths can be updated on command line using the following parameters:
- `--fasta` - Path to genome fasta file
- `--bowtie2Index` - Path to HiSAT2 index
- `--bwaIndex` - Path to TopHat2 index
- `--gtf` - Path to GTF file
- `--geneBed` - Path to gene file
- `--effGenomeSize` - Effective Genome Size
- `--blacklist` - Path to black list region file (bed format)
- `--saveAlignedIntermediates` - Save the BAM files from the Aligment step  - not done by default


### `--genomeAnnotationPath`

The `--genomeAnnotationPath` define where all annotations are stored. This path can be defined on the command line or set up in the different configuration file during the pipeline installation.
See `conf/installation.md` for details.

## Alignment

### `--aligner`

Specify which tool must be used for reads alignment. The expected values are `bwa-mem` or `bowtie2` (default: `bwa-mem`).

```bash
--aligner 'bowtie2' 
```

## Filtering

### `--mapq`

Filter all reads in the alignment files with a mapping quality lower than this threshold.
By default, no mapq filtering is performed.

```bash
--mapq 20
```

### `--keepDups`

By default, duplicates are removed from the aligned data for downstream analysis.
Use this option to keep the duplicates.

```bash
--keepDuplicates
```

### `--keepMito`

It is usually recommanded to remove reads aligned on mitochondrial chromosome which can represent 
a high fraction of reads for some cell types.
Use this option if you do not want to apply this filter.

```bash
--keepMito
```

### `--keepSingleton`

Use this option to keep singleton reads in the analysis which are usually removed by default.

```bash
--keepSingleton
```

## Peak calling

### `--caller`

Specify which tool(s) to use to call peaks in the ATAC-seq data. Available tools are `macs2` and `genrich`.
Both tools can be specified (comma separated).

```bash
--caller 'macs2,genrich' 
```

### `--tn5sites`

Use this option to focus the analysis on the transposase insertion site, which are usually centered around 5' end
of R1/R2 reads. Activating this option has an impact on `bigwig` files and peak calling parameters.
This parameter is usually advised for motif analysis such as footprinting.

When the option is turned off (the default), reads are usually extended in order to work at the fragment level.

```bash
--tn5sites
```

### `--extsize`

Extension parameter for macs2 peak calling when `--tn5sites` is specified.
The macs2 shift parameter is automatically set as extsize/2.

```
--extisze 150
```

## Profiles

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile docker` - the order of arguments is important!

The following `-profile` are available. If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

  - `test`
  
  A profile with a complete configuration for automated testing. It includes links to test data so needs no other parameters.

  - `conda`
  
  Build a new conda environment before running the pipeline. Use the option `--condaCacheDir` to change the default conda cache directory.
  
  - `multiconda`
  
  Build a new conda environment for each process before running the pipeline. Use the option ``--condaCacheDir` to change the default conda cache directory.

  - `path`
  
  Use a global path for all tools. Use the option `--globalPath` to define the path the use.
  
  - `multipath`
  
  Use the paths defined in configuration for each tool.
  
  - `docker`
  
  Use the Docker images for each process.
  
  - `singularity`
  
  Use the Singularity images for each process. Use the option `--singularityImagePath` to specify where the images are available.
  
  - `cluster`
  
  Submit the jobs on the cluster instead of running them locally.
												

## Other command line parameters

### `--skip*`

The pipeline is made with a few *skip* options that allow to skip optional steps in the workflow.
The following options can be used:
- `--skipFastqc` - Skip FastQC
- `--skipPreseq` - Skips preseq QC
- `--skipShift` - Skips reads shifting for Tn5 correction (+4/-5bp)
- `--skipDeepTools` - Skips deeptools QC
- `--skipPeakCalling` - Skips peak calling
- `--skipPeakAnno` - Skips peak annotation
- `--skipMultiqc` - Skip MultiQC
				
### `--metadata`

Specify a two-columns (tab-delimited) metadata file to diplay in the final Multiqc report.

### `--outDir`

The output directory where the results will be saved.

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--maxMemory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--maxMemory '8.GB'`

### `--maxTime`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--maxTime '2.h'`

### `--maxCpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--maxCpus 1`

### `--multiqcConfig`

Specify a path to a custom MultiQC configuration file.

## Job resources

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time (see the `conf/base.conf` file). 
For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.
