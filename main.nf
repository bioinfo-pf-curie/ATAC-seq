#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2020
This software is a computer program whose purpose is to analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data.
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

This script is based on the nf-core guidelines. See https://nf-co.re/ for more information
*/


/*
========================================================================================
            ATAC-seq
========================================================================================
ATAC-seq Analysis Pipeline.
#### Homepage / Documentation
https://gitlab.curie.fr/data-analysis/atac-seq
----------------------------------------------------------------------------------------
*/

def helpMessage() {
  if ("${workflow.manifest.version}" =~ /dev/ ){
  devMess = file("$baseDir/assets/dev_message.txt")
  log.info devMess.text
  }

  log.info"""

  ATAC-seq v${workflow.manifest.version}
  ======================================================================

  Usage:

  nextflow run main.nf --reads '*_R{1,2}.fastq.gz' -profile conda --genomeAnnotationPath '/data/annotations/pipelines' --genome 'hg19' 
  nextflow run main.nf --samplePlan 'sample_plan.csv' --design 'design.csv' -profile conda --genomeAnnotationPath '/data/annotations/pipelines' --genome 'hg19'

  Mandatory arguments:
  --reads [file]                     Path to input data (must be surrounded with quotes)
  --samplePlan [file]                Path to sample plan file if '--reads' is not specified
  --genome [str]                     Name of genome reference. See the `--genomeAnnotationPath` to defined the annotations path.
  -profile [str]                     Configuration profile to use. Can use multiple (comma separated)

  Inputs:
  --singleEnd [bool]                 Specifies that the input is single end reads
  --fragmentSize [int]               Estimated fragment length used to extend single-end reads. Default: 0

  References           If not specified in the configuration file or you wish to overwrite any of the references given by the --genome field
  --fasta [file]                     Path to Fasta reference

  Alignment:
  --aligner [str]                    Alignment tool to use ['bwa-mem', 'bowtie2']. Default: 'bwa-mem'
  --saveAlignedIntermediates [bool]  Save all intermediates mapping files. Default: false  
  --bwaIndex [file]                  Index for Bwa-mem aligner
  --bowtie2Index [file]              Index for Bowtie2 aligner
  --bwaOpts [str]                    Modify the Bwa-mem mapping parameters
  --bowtie2Opts [str]                Modify the Bowtie2 mapping parameters

  Trimming:
  --trimNextseq [int]            Instructs Trim Galore to apply the --nextseq=X option, to trim based on quality after removing poly-G tails (Default: 0)
  --skipTrimming [bool]          Skip the adapter trimming step (Default: false)
  --saveTrimmed [bool]           Save the trimmed FastQ files in the results directory (Default: false)

  Filtering:
  --mapq [int]                       Minimum mapping quality to consider. Default: 0
  --keepDups [bool]                  Do not remove duplicates afer marking. Default: false
  --keepSingleton [bool]             Keep unpaired reads. Default: false
  --keepMito [bool]                  Do not filter reads from mitochrondrial chromosomal. Default: false
  --blacklist [file]                 Path to black list regions (.bed).
  --ignoreBlacklist [bool]           Do not filter blacklisted regions. Default: false

  Calling:
  --caller [str]                     Peak caller to use ['macs2','genrich']. Several tools can be specified (comma separated). Default: 'macs2'
  --tn5sites [bool]                  Focus the analysis on Tn5 insertion sites (ie. work at the reads level and not at the fragment one). Default: false
  --extsize [int]                    Value to use for extsize parameter during Macs calling. Shift parameter will be set up as extsize/2. Default: 73

  Annotation:          If not specified in the configuration file or you wish to overwrite any of the references given by the --genome field
  --genomeAnnotationPath             Path to genome annotations.
  --geneBed [file]                   BED annotation file with gene coordinate.
  --gtf [file]                       GTF annotation file. Used in HOMER peak annotation
  --effGenomeSize [int]              Effective Genome size

  Skip options:        All are false by default
  --skipFastqc [bool]                Skips fastQC. Default: false
  --skipPreseq [bool]                Skips preseq QC. Default: false
  --skipShift [bool]                 Skips reads shifting for Tn5 correction (+4/-5bp). Default: false
  --skipDeepTools [bool]             Skips deeptools QC. Default: false
  --skipPeakCalling [bool]           Skips peak calling. Default: false
  --skipPeakAnno                     Skips peak annotation. Default: false
  --skipMultiQC [bool]               Skips MultiQC step. Default: false

  Other options:
  --metadata [file]                  Path to metadata file for MultiQC report
  --outDir [file]                    The output directory where the results will be saved
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

  """.stripIndent()
}

// Show help messsage
if (params.help){
  helpMessage()
  exit 0
}

if (params.aligner != 'bwa-mem' && params.aligner != 'star' && params.aligner != 'bowtie2' ) {
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star', 'bowtie2' or 'bwa-mem'"
}

/*********************
 * Fasta file
 */

// Configurable reference genomes
genomeRef = params.genome

// Genome Fasta file
params.fasta = genomeRef ? params.genomes[ genomeRef ].fasta ?: false : false
if ( params.fasta ){
  Channel
    .fromPath(params.fasta, checkIfExists: true)
    .set{chFastaHomer}
}
else{
  exit 1, "Fasta file not found: ${params.fasta}"
}

// Chromosome size file
params.chrsize = genomeRef ? params.genomes[ genomeRef ].chrsize ?: false : false
if ( params.chrsize ){
  Channel
    .fromPath(params.chrsize, checkIfExists: true)
    .set{chChromSize}
}
else{
  exit 1, "Chromosome size file not found: ${params.chrsize}"
}

params.mitoName = genomeRef ? params.genomes[ genomeRef ].mitoName ?: false : false


/********************
 * Bwa-mem Index
 */

params.bwaIndex = genomeRef ? params.genomes[ genomeRef ].bwaIndex ?: false : false
if (params.bwaIndex){
  lastPath = params.bwaIndex.lastIndexOf(File.separator)
  bwaDir =  params.bwaIndex.substring(0,lastPath+1)
  bwaBase = params.bwaIndex.substring(lastPath+1)
  Channel
    .fromPath(bwaDir, checkIfExists: true)
    .ifEmpty {exit 1, "BWA index file not found: ${params.bwaIndex}"}
    .combine( [ bwaBase ] )
    .dump(tag :'bwaindexch')
    .set { chBwaIndex }
} else {
  exit 1, "BWA index file not found: ${params.bwaIndex}"
}


/*********************
 * Bowtie2 indexes
 */

params.bt2Index = genomeRef ? params.genomes[ genomeRef ].bowtie2Index ?: false : false
if (params.bt2Index){
  lastPath = params.bt2Index.lastIndexOf(File.separator)
  bt2Dir =  params.bt2Index.substring(0,lastPath+1)
  bt2Base = params.bt2Index.substring(lastPath+1)
  Channel
    .fromPath(bt2Dir, checkIfExists: true)
    .ifEmpty {exit 1, "Bowtie2 index file not found: ${params.bt2Index}"}
    .combine( [ bt2Base ] ) 
    .set { chBt2Index }
} else {
  exit 1, "Bowtie2 index file not found: ${params.bt2Index}"
}

/*********************
 * Annotations
 */

params.gtf = genomeRef ? params.genomes[ genomeRef ].gtf ?: false : false
if (params.gtf) {
  Channel
    .fromPath(params.gtf, checkIfExists: true)
    .set{chGtfHomer}
}
else {
  exit 1, "GTF annotation file not specified!"
}

params.geneBed = genomeRef ? params.genomes[ genomeRef ].geneBed ?: false : false
if (params.geneBed) {
  Channel
    .fromPath(params.geneBed, checkIfExists: true)
    .ifEmpty {exit 1, "BED file ${geneBed} not found"}
    .into{chGeneBed; chGeneBedDeeptools; chGenePrepareAnnot}
}else{
  chGeneBed = Channel.empty()
  chGeneBedDeeptools = Channel.empty()
  chGenePrepareAnnot = Channel.empty()
}

params.blacklist = genomeRef ? params.genomes[ genomeRef ].blacklist ?: false : false
if (params.blacklist && !params.ignoreBlacklist) { 
  Channel
    .fromPath(params.blacklist, checkIfExists: true)
    .into {chBlacklistBigWig; chBlacklistCorrelation} 
}else{
  chBlacklistBigWig = Channel.empty()
  chBlacklistCorrelation = Channel.empty()
}


/***********************
 * Header and conf
 */

//Peak Calling
params.effGenomeSize = genomeRef ? params.genomes[ genomeRef ].effGenomeSize ?: false : false
if (!params.effGenomeSize) {
  log.warn "=================================================================\n" +
            "  WARNING! Effective Genome Size is not defined.\n" +
            "  Peak calling and annotation will be skipped.\n" +
            "  Please specify value for '--effGenomeSize' to run these steps.\n" +
            "================================================================"
}

Channel
  .fromPath("$baseDir/assets/peak_count_header.txt")
  .into { chPeakCountHeaderMacs; chPeakCountHeaderGenrich }

Channel
  .fromPath("$baseDir/assets/frip_score_header.txt")
  .into { chFripScoreHeaderMacs; chFripScoreHeaderGenrich }


Channel
  .fromPath("$baseDir/assets/peak_annotation_header.txt")
  .set{ chPeakAnnotationHeader }

//Stage config files
Channel
  .fromPath(params.multiqcConfig, checkIfExists: true)
  .set{chMultiqcConfig}
chOutputDocs = file("$baseDir/docs/output.md", checkIfExists: true)
chOutputDocsImages = file("$baseDir/docs/images/", checkIfExists: true)

//Has the run name been specified by the user?
//This has the bonus effect of catching both -name and --name
customRunName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  customRunName = workflow.runName
}

//Header log info
if ("${workflow.manifest.version}" =~ /dev/ ){
  devMess = file("$baseDir/assets/dev_message.txt")
  log.info devMess.text
}

log.info """=======================================================

ATAC-seq v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'ATAC-seq'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = customRunName ?: workflow.runName
if (params.samplePlan) {
  summary['SamplePlan'] = params.samplePlan
} else{
  summary['Reads']      = params.reads
}
summary['Annotation']   = params.genomeAnnotationPath
summary['Fasta Ref']    = params.fasta
summary['GTF']          = params.gtf
summary['Genes']        = params.geneBed
if (params.blacklist && !params.ignoreBlacklist)  summary['Blacklist '] = params.blacklist
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
if (params.singleEnd)  summary['Fragment Size '] = params.fragmentSize
summary['Aligner'] = params.aligner
summary['Mode']         = params.tn5sites ? 'Tn5-sites' : 'Fragment'
summary['Keep Duplicates'] = params.keepDups ? 'Yes':'No'
summary['Keep Singleton'] = params.keepSingleton ? 'Yes':'No'
summary['Keep Mito'] = params.keepMito ? 'Yes':'No'
summary['Min MapQ'] = params.mapq
summary['Peak Calling']= params.caller
summary['Max Memory']   = params.maxMemory
summary['Max CPUs']     = params.maxCpus
summary['Max Time']     = params.maxTime
summary['Output dir']   = params.outDir
summary['Working dir']  = workflow.workDir
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outDir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile

if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

/*
 * CHANNELS
 */

if ( params.metadata ){
  Channel
    .fromPath( params.metadata )
    .ifEmpty { exit 1, "Metadata file not found: ${params.metadata}" }
    .set { chMetadata }
}

/*
 * Create a channel for input read files
 */

if(params.samplePlan){
   if(params.singleEnd && !params.inputBam){
      Channel
         .from(file("${params.samplePlan}"))
         .splitCsv(header: false)
         .map{ row -> [ row[0], [file(row[2])]] }
         //.join(chDesignCsv2)
         //.map { row -> [ row[2] + '_' + row[4], [file(row[1][0])], row[0], row[3], row[4] ] }
         .into { rawReadsFastqc; rawReads; planMultiQC }

   }else if (!params.singleEnd && !params.inputBam){
      Channel
         .from(file("${params.samplePlan}"))
         .splitCsv(header: false)
         .map{ row -> [ row[0], [file(row[2]), file(row[3])]] }
         //.join(chDesignCsv2)
         //.map { row -> [ row[2] + '_' + row[4], [file(row[1][0]),file(row[1][1])], row[0], row[3], row[4] ] }
         .into { rawReadsFastqc; rawReads; planMultiQC }

   }else{
      Channel
         .from(file("${params.samplePlan}"))
         .splitCsv(header: false)
         .map{ row -> [ row[0], [file(row[2])]]}
         .dump(tag: 'alignedreads')
         //.join(chDesignCsv2)
         //.map { row -> [ row[2] + '_' + row[4], [file(row[1][0])], row[0], row[3], row[4] ] }
         .into { chAlignReads; planMultiQC }
       params.reads=false

       Channel.empty().into{ rawReadsFastqc; rawReads }
  }
} else if(params.readPaths){
    if(params.singleEnd){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            //.join(chDesignCsv2)
            //.map { row -> [ row[2] + '_' + row[4], [file(row[1][0])], row[0], row[3], row[4] ] }
            .into { rawReadsFastqc; rawReads; planMultiQC }

    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            //.join(chDesignCsv2)
            //.map { row -> [ row[2] + '_' + row[4], [file(row[1][0]),file(row[1][1])], row[0], row[3], row[4] ] }
            .into { rawReadsFastqc; rawReads; planMultiQC }
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        //.join(chDesignCsv2)
        //.map { row -> [ row[2] + '_' + row[4], [file(row[1][0]),file(row[1][1])], row[0], row[3], row[4] ] }
        .into { rawReadsFastqc; rawReads; planMultiQC }
}


/**************************
 * Make sample plan if not available
 */

if (params.samplePlan){
  Channel
    .fromPath(params.samplePlan)
    .into { chSplan; chSplanCheck }
}else if(params.readPaths){
  if (params.singleEnd){
    Channel
       .from(params.readPaths)
       .collectFile() {
         item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
        }
       .into { chSplan; chSplanCheck }
  }else{
     Channel
       .from(params.readPaths)
       .collectFile() {
         item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
        }
       .into { chSplan; chSplanCheck }
  }
} else {
  if (params.singleEnd){
    Channel
       .fromFilePairs( params.reads, size: 1 )
       .collectFile() {
          item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n']
       }     
       .into { chSplan; chSplanCheck }
  }else{
    Channel
       .fromFilePairs( params.reads, size: 2 )
       .collectFile() {
          item -> ["sample_plan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
       }     
       .into { chSplan; chSplanCheck }
   }
}


/*********************************************************
/*
 * FastQC
 */

process fastQC{
  tag "${prefix}"
  label 'fastqc'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outDir}/fastqc", mode: 'copy'

  when:
  !params.skipFastqc && !params.inputBam

  input:
  set val(prefix), file(reads) from rawReadsFastqc

  output:
  file "*_fastqc.{zip,html}" into chFastqcMqc
  file("v_fastqc.txt") into chFastqcVersion

  script:
  """
  fastqc --version &> v_fastqc.txt
  fastqc -q $reads -t ${task.cpus}
  """
}

/********************************************************
/*
 * Trim Galore 
 */

if (params.skipTrimming) {
    rawReads  
        .into {trimmedReadsBWA ; trimmedReadsBt2 }
} else {
process trimgalore {
  tag "${prefix}"
  label 'trimgalore'
  label 'highCpu'
  label 'highMem'
  publishDir "${params.outDir}/trimgalore", mode: 'copy'

  when:
  !params.inputBam && !params.skipTrimming

  input:
  set val(prefix), file(reads) from rawReads

  output:
  tuple val(prefix), path('*.fq.gz') into trimmedReadsBWA,trimmedReadsBt2

  script:
    nextSeq = params.trimNextseq > 0 ? "--nextseq ${params.trimNextseq}" : ''

  // Added soft-links to original fastqs for consistent naming in MultiQC
    if (params.singleEnd) {
      """
      [ ! -f  ${prefix}.fastq.gz ] && ln -s ${reads} ${prefix}.fastq.gz
      trim_galore --cores ${task.cpus} --fastqc --gzip $nextSeq ${prefix}.fastq.gz
      """
    } else {
      """
      [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
      [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
      trim_galore --cores ${task.cpus} --paired --fastqc --gzip $nextSeq ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz
      """
    }
  }
}



/*
 * Alignment on reference genome
 */

/* BWA-MEM */
process bwaMem{
  tag "${sample} on ${genomeBase}"
  label 'bwa'
  label 'highCpu'
  label 'extraMem'
  publishDir "${params.outDir}/mapping", mode: 'copy',
             saveAs: {filename -> 
	     if (filename.indexOf(".log") > 0) "logs/$filename" 
	     else if (params.saveAlignedIntermediates) filename}

  when:
  params.aligner == "bwa-mem" && !params.inputBam

  input:
  set val(sample), file(reads), file(index), val(genomeBase) from trimmedReadsBWA.combine(chBwaIndex)

  output:
  set val(sample), file("*.bam") into chAlignReadsBwa
  file("*.log") into chBwaMqc
  file("v_bwa.txt") into chBwaVersion

  script:
  """
  echo \$(bwa 2>&1) &> v_bwa.txt
  bwa mem -t ${task.cpus} \\
           ${index}/${genomeBase} \\
          ${params.bwaOpts} \\
          $reads | samtools view -bS - > ${sample}.bam
  ${baseDir}/bin/getBWAstats.sh -i ${sample}.bam > ${sample}_bwa.log
  """
}

/* BOWTIE2 */
process bowtie2{
  tag "${sample} on ${genomeBase}"
  label 'bowtie2'
  label 'highCpu'
  label 'medMem'
  publishDir "${params.outDir}/mapping", mode: 'copy',
              saveAs: {filename ->
	      if (filename.indexOf(".log") > 0) "logs/$filename"  
	      else if (params.saveAlignedIntermediates) filename}
  when:
  params.aligner == "bowtie2" && !params.inputBam

  input:
  set val(sample), file(reads), file(index), val(genomeBase) from trimmedReadsBt2.combine(chBt2Index)

  output:
  set val(sample), file("*.bam") into chAlignReadsBowtie2
  file("*.log") into chBowtie2Mqc
  file("v_bowtie2.txt") into chBowtie2Version

  script:
  readCommand = params.singleEnd ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
  """
  bowtie2 --version &> v_bowtie2.txt
  bowtie2 -p ${task.cpus} \\
          ${params.bowtie2Opts} \\
           -x ${index}/${genomeBase} \\
          $readCommand > ${sample}.bam 2> ${sample}_bowtie2.log
  """
}

if(!params.inputBam){
if (params.aligner == "bowtie2"){
  chAlignReads = chAlignReadsBowtie2
  chMappingMqc = chBowtie2Mqc
} else if (params.aligner == "bwa-mem"){
  chAlignReads = chAlignReadsBwa
  chMappingMqc = chBwaMqc
} else if (params.aligner == "star"){
  chAlignReads = chAlignReadsStar
  chMappingMqc = chStarMqc
}
}

if (params.inputBam){
  chFastqcMqc = Channel.empty()
  chMappingMqc = Channel.empty()
}

/*
 * Sorting BAM files
 */

process bamSort{
  tag "${prefix}"
  label 'samtools'
  label 'medCpu'
  label 'medMem'
  publishDir path: "${params.outDir}/mapping", mode: 'copy',
    saveAs: {filename ->
             if ( filename.endsWith("stats") && params.saveAlignedIntermediates ) "stats/$filename"
             else if ( (filename.endsWith(".bam") || (filename.endsWith(".bam.bai"))) && params.saveAlignedIntermediates ) filename
             else null
            }

  input:
  set val(prefix), file(unsortedBam) from chAlignReads

  output:
  set val(prefix), file('*sorted.{bam,bam.bai}') into chSortBams
  file("*stats") into chStats
  file("*mqc") into chStatsMqc
  file("v_samtools.txt") into chSamtoolsVersionBamSort

  script:
  """
  samtools --version &> v_samtools.txt
  samtools sort $unsortedBam -@ ${task.cpus} -T ${prefix} -o ${prefix}_sorted.bam
  samtools index ${prefix}_sorted.bam
  samtools flagstat ${prefix}_sorted.bam > ${prefix}_sorted.flagstats
  samtools idxstats ${prefix}_sorted.bam > ${prefix}_sorted.idxstats
  samtools stats ${prefix}_sorted.bam > ${prefix}_sorted.stats

  aligned="\$(samtools view -@ $task.cpus -F 0x100 -F 0x4 -F 0x800 -c ${prefix}_sorted.bam)"
  hqbam="\$(samtools view -@ $task.cpus -F 0x100 -F 0x800 -F 0x4 -q 20 -c ${prefix}_sorted.bam)"
  lqbam="\$((\$aligned - \$hqbam))"
  echo -e "Mapped,\${aligned}" > ${prefix}_mappingstats.mqc
  echo -e "HighQual,\${hqbam}" >> ${prefix}_mappingstats.mqc
  echo -e "LowQual,\${lqbam}" >> ${prefix}_mappingstats.mqc
  """
}

/*
 * Marking duplicates
 */

process markDuplicates{
  tag "${prefix}"
  label 'picard'
  label 'medCpu'
  label 'medMem'
  publishDir path: "${params.outDir}/mapping", mode: 'copy',
    saveAs: {filename ->
             if (!filename.endsWith(".bam") && !filename.endsWith(".bam.bai") && params.saveAlignedIntermediates ) "stats/$filename"
             else if ( (filename.endsWith(".bam") || (filename.endsWith(".bam.bai"))) && params.saveAlignedIntermediates ) filename
             else null}

  input:
  set val(prefix), file(sortedBams) from chSortBams

  output:
  set val(prefix), file("*marked.{bam,bam.bai}") into chMarkedBams, chMarkedBamsFilt, chMarkedPreseq
  set val(prefix), file("*marked.flagstat") into chMarkedFlagstat
  file "*marked.{idxstats,stats}" into chMarkedStats
  file "*metrics.txt" into chMarkedPicstats
  file("v_picard.txt") into chPicardVersion

  script:
  """
  echo \$(picard MarkDuplicates --version 2>&1) &> v_picard.txt
  picard -Xmx4g MarkDuplicates \\
    INPUT=${sortedBams[0]} \\
    OUTPUT=${prefix}_marked.bam \\
    ASSUME_SORTED=true \\
    REMOVE_DUPLICATES=false \\
    METRICS_FILE=${prefix}.MarkDuplicates.metrics.txt \\
    VALIDATION_STRINGENCY=LENIENT \\
    TMP_DIR=tmpdir
  samtools index ${prefix}_marked.bam
  samtools idxstats ${prefix}_marked.bam > ${prefix}_marked.idxstats
  samtools flagstat ${prefix}_marked.bam > ${prefix}_marked.flagstat
  samtools stats ${prefix}_marked.bam > ${prefix}_marked.stats
  """
}

/*
 * Preseq (before alignment filtering)
 */

process preseq {
  tag "${prefix}"
  label 'preseq'
  label 'lowCpu'
  label 'medMem'
  publishDir "${params.outDir}/preseq", mode: 'copy'

  when:
  !params.skipPreseq

  input:
  set val(prefix), file(bam) from chMarkedPreseq

  output:
  file "*.ccurve.txt" into chPreseqStats
  file("v_preseq.txt") into chPreseqVersion

  script:
  defectMode = params.preseqDefect ? '-D' : ''
  """
  preseq &> v_preseq.txt
  preseq lc_extrap -v $defectMode -output ${prefix}.ccurve.txt -bam ${bam[0]} -e  300e+06
  """
}


/*
 * BAM Filtering
 */

process bamFiltering {
  tag "${prefix}"
  label 'samtools'
  label 'lowCpu'
  label 'medMem'
  publishDir path: "${params.outDir}/mapping", mode: 'copy',
    saveAs: {filename ->
             if (!filename.endsWith(".bam") && (!filename.endsWith(".bam.bai"))) "stats/$filename"
             else if ((filename.endsWith("_filtered.bam") || (filename.endsWith("_filtered.bam.bai"))) && (params.saveAlignedIntermediates || params.skipShift)) filename
             else null}

  input:
  set val(prefix), file(markedBam) from chMarkedBamsFilt
  file bed from chGeneBed.collect()

  output:
  set val(prefix), file("*filtered.{bam,bam.bai}") into chFilteredBams,chFilteredBamsNoShift
  set val(prefix), file("*filtered.flagstat") into chFilteredFlagstat
  file "*raw.idxstats" into chRawStats
  file "*filtered.{idxstats,stats}" into chFilteredStats
  file("v_samtools.txt") into chSamtoolsVersionBamFiltering

  script:
  filterParams = params.singleEnd ? "-F 0x004" : params.keepSingleton ? "-F 0x004" : "-F 0x004 -F 0x0008 -f 0x001"
  dupParams = params.keepDups ? "" : "-F 0x0400"
  mapqParams = params.mapq > 0 ? "-q ${params.mapq}" : ""
  mitoParams = params.keepMito ? "" : params.mitoName ? "grep -v ${params.mitoName} |" : ""
  nameSortBam = params.singleEnd ? "" : "samtools sort -n -@ $task.cpus -o ${prefix}.bam -T $prefix ${prefix}_filtered.bam"
  """
  samtools --version &> v_samtools.txt
  samtools idxstats ${markedBam[0]}  >  ${prefix}_raw.idxstats
  cat ${prefix}_raw.idxstats | cut -f 1 | ${mitoParams} xargs samtools view \\
    $filterParams \\
    $dupParams \\
    $mapqParams \\
    -b ${markedBam[0]} > ${prefix}_filtered.bam
  samtools index ${prefix}_filtered.bam
  samtools flagstat ${prefix}_filtered.bam > ${prefix}_filtered.flagstat
  samtools idxstats ${prefix}_filtered.bam > ${prefix}_filtered.idxstats
  samtools stats ${prefix}_filtered.bam > ${prefix}_filtered.stats
  $nameSortBam
  """
}


/*
 * Shifting reads
 */

process readsShifting {
  tag "${prefix}"
  label 'deeptools'
  label 'medCpu'
  label 'medMem'
  publishDir path: "${params.outDir}/mapping", mode: 'copy'

  when: !params.skipShift

  input:
  set val(prefix), file(filteredBam) from chFilteredBams 

  output:
  set val(prefix), file("*filtered_shifted.{bam,bam.bai}") into chShiftBams

  script:
  """
  alignmentSieve -b ${filteredBam[0]} -p ${task.cpus} --ATACshift -o ${prefix}_filtered_shifted_tmp.bam
  samtools sort  ${prefix}_filtered_shifted_tmp.bam -@ ${task.cpus} -o ${prefix}_filtered_shifted.bam 
  samtools index ${prefix}_filtered_shifted.bam
  """
}

if (! params.skipShift){
  chShiftBams
    .dump(tag : 'fbams')
    .into{ chBamsFragSize;
           chBamsMacs; chBamsGenrich;
           chBamsBigWig;
           chBamDTCor ; chBaiDTCor; chSampleDTCor ;
           chBamDTFingerprint ; chBaiDTFingerprint ; chSampleDTFingerprint }
}else{
  chFilteredBamsNoShift
    .dump(tag : 'fbams')
    .into{ chBamsFragSize;
           chBamsMacs; chBamsGenrich;
           chBamsBigWig;
           chBamDTCor ; chBaiDTCor; chSampleDTCor ;
           chBamDTFingerprint ; chBaiDTFingerprint ; chSampleDTFingerprint }
}

/*
 * Get fragment sizes
 */

process getFragmentSize {
  tag "${prefix}"
  label 'picard'
  label 'lowCpu'
  label 'medMem'

  publishDir path: "${params.outDir}/fragSize/", mode: "copy"

  when:
  !params.singleEnd 

  input:
  set val(prefix), file(filteredBam) from chBamsFragSize

  output:
  file("*.{pdf,txt}") into chFragmentsSize

  script:
  """
  picard CollectInsertSizeMetrics \
      I=${filteredBam[0]} \
      O=${prefix}_insert_size_metrics.txt \
      H=${prefix}_insert_size_histogram.pdf \
      VALIDATION_STRINGENCY=LENIENT \
      M=0.5
  """
}

/* 
 * BigWig Tracks
 */

process bigWig {
  tag "${prefix}"
  label 'deeptools'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outDir}/bigWig", mode: "copy",
              saveAs: {filename -> if (filename.endsWith("norm.bigwig")) filename
                       else null}

  input:
  set val(prefix), file(filteredBams) from chBamsBigWig
  file(BLbed) from chBlacklistBigWig.collect().ifEmpty([])

  output:
  set val(prefix), file('*norm.bigwig') into chBigWig
  set val("${prefix}_NFR"), file('*NFR.bigwig') into chBigWigNFR
  set val("${prefix}_NBR"), file('*NBR.bigwig') into chBigWigNBR
  file("v_deeptools.txt") into chDeeptoolsVersion

  script:
  extend = params.tn5sites ? "" : params.singleEnd && params.fragmentSize > 0 ? "--extendReads ${params.fragmentSize}" : "--extendReads"
  blacklistParams = params.blacklist && !params.ignoreBlacklist ? "--blackListFileName ${BLbed}" : ""
  effGsize = params.effGenomeSize ? "--effectiveGenomeSize ${params.effGenomeSize}" : ""
  """
  bamCoverage --version &> v_deeptools.txt
  nbreads=\$(samtools view -@ $task.cpus -F 0x100 -F 0x4 -F 0x800 -c ${filteredBams[0]})
  sf=\$(echo "1000000 \$nbreads" | awk '{print \$1/\$2}')

  bamCoverage -b ${filteredBams[0]} \\
              -o ${prefix}_norm.bigwig \\
              -p ${task.cpus} \\
              ${blacklistParams} \\
              ${effGsize} \\
	      ${extend} \\
              --scaleFactor \$sf

  bamCoverage -b ${filteredBams[0]} \\
              -o ${prefix}_NFR.bigwig \\
              -p ${task.cpus} \\
              ${blacklistParams} \\
              ${effGsize} \\
              --maxFragmentLength 100 \\
              --scaleFactor \$sf

  bamCoverage -b ${filteredBams[0]} \\
              -o ${prefix}_NBR.bigwig \\
              -p ${task.cpus} \\
              ${blacklistParams} \\
              ${effGsize} \\
              --minFragmentLength 180 \\
              --maxFragmentLength 473 \\
              --scaleFactor \$sf
  """
}

/*
 * DeepTools QC
 */

process deepToolsComputeMatrix{
  tag "${prefix}"
  label 'deeptools'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outDir}/deepTools/computeMatrix", mode: "copy"

  input:
  set val(prefix), file(bigwig) from chBigWigNFR.concat(chBigWigNBR)
  file geneBed from chGeneBedDeeptools.collect()

  output:
  file("*.{gz,pdf,png}") optional true into chDeeptoolsSingle
  file("*corrected.tab") optional true into chDeeptoolsSingleMqc

  when:
  !params.skipDeepTools 

  script:
  effGsize = params.effGenomeSize ? "--effectiveGenomeSize ${params.effGenomeSize}" : ""
  """
  computeMatrix reference-point \\
                --referencePoint TSS \\
                -R $geneBed \\
                -S ${bigwig} \\
                -o ${prefix}_matrix.mat.gz \\
                --outFileNameMatrix ${prefix}.computeMatrix.vals.mat.gz \\
                -b 1000 -a 1000 --skipZeros \\
                -p ${task.cpus} --sortRegions no

  ## if matrix is not empty
  if [ `zcat ${prefix}_matrix.mat.gz | wc -l` -ge "2" ]
  then
    plotProfile -m ${prefix}_matrix.mat.gz \\
              -o ${prefix}_profile.pdf \\
              --outFileNameData ${prefix}.plotProfile.tab
  
    sed -e 's/.0\t/\t/g' -e 's/tick/TSS/' ${prefix}.plotProfile.tab | sed -e 's@.0\$@@g' > ${prefix}.plotProfile_corrected.tab

    plotHeatmap -m ${prefix}_matrix.mat.gz -out ${prefix}_plotHeatmap.png
  fi
  """
}

process deepToolsFingerprint{
  label 'deeptools'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outDir}/deepTools/fingerprintQC", mode: "copy"

  when:
  !params.skipDeepTools

  input:
  file(allBams) from chBamDTFingerprint.map{it[1][0]}.collect()
  file(allBai) from chBaiDTFingerprint.map{it[1][1]}.collect()
  val (allPrefix) from chSampleDTFingerprint.map{it[0]}.collect()
 
  output:
  file "bams_fingerprint.pdf" into chDeeptoolsFingerprint
  file "plotFingerprint*" into chDeeptoolsFingerprintMqc 
 
  script:
  extend = params.tn5sites ? "" : params.singleEnd && params.fragmentSize > 0 ? "--extendReads ${params.fragmentSize}" : "--extendReads"
  allPrefix = allPrefix.toString().replace("[","")
  allPrefix = allPrefix.replace(","," ") 
  allPrefix = allPrefix.replace("]","")
  """
  plotFingerprint -b $allBams \\
                  -plot bams_fingerprint.pdf \\
                  -p ${task.cpus} \\
                  -l $allPrefix \\
		  ${extend} \\
                  --skipZeros \\
                  --outRawCounts plotFingerprint.raw.txt \\
                  --outQualityMetrics plotFingerprint.qmetrics.txt
  """
}


/***********************
 * Peak calling 
 */

(chFilteredMacsFlagstat, chFilteredGenrichFlagstat) = chFilteredFlagstat.into(2)

/*
 * MACS2
 */

process macs2 {
  tag "${prefix}"
  label 'macs2'
  label 'medCpu'
  label 'medMem'
  publishDir path: "${params.outDir}/peakCalling/macs2", mode: 'copy',
    saveAs: { filename -> if (filename.endsWith(".tsv")) "stats/$filename"
                          else filename }
 
  when:
  !params.skipPeakCalling && params.caller =~ "macs2"

  input:
  set val(prefix), file(bam), file(sampleFlagstat) from chBamsMacs.join(chFilteredMacsFlagstat)
  file peakCountHeader from chPeakCountHeaderMacs.collect()
  file fripScoreHeader from chFripScoreHeaderMacs.collect()

  output:
  file("*macs2*.{xls,bed}") into chMacsOutput
  set val(prefix), file("*macs2_peaks.narrowPeak"), val("macs2") into chMacsPeaks,chMacsPeaksbb
  set val(prefix), file("*macs2_summits*") into chMacs2Summits
  file "*_mqc.tsv" into chMacsMqc
  file("v_macs2.txt") into chMacsVersion
  file("v_bedtools.txt") into chMacsBedtoolsVersion

  script:
  preproc = params.tn5sites ? "bamToBed -i ${bam[0]} > ${bam[0].baseName}.bed; shift=\$(expr ${params.extsize} / 2 )" : ""
  inputs = params.tn5sites ? "-t ${bam[0].baseName}.bed -f BED" : params.singleEnd ? "-t ${bam[0]}-f BAM" : "-t ${bam[0]} -f BAMPE"
  opts = params.tn5sites ? "--keep-dup all --shift -\${shift} --extsize ${params.extsize} --nomodel -q 0.01" : "--keep-dup all --nomodel -q 0.01"
  """
  echo \$(macs2 --version 2>&1) &> v_macs2.txt
  echo \$(bedtools --version 2>&1) &> v_bedtools.txt
  ${preproc}
  macs2 callpeak \\
    ${inputs} \\
    -g $params.effGenomeSize \\
    -n ${prefix}_macs2 \\
    --SPMR --trackline --bdg \\
    ${opts}

  cat ${prefix}_macs2_peaks.narrowPeak | tail -n +2 | wc -l | awk -v OFS='\t' '{ print "${prefix}", \$1 }' | cat $peakCountHeader - > ${prefix}_macs2_peaks.count_mqc.tsv
  READS_IN_PEAKS=\$(intersectBed -a ${bam[0]} -b ${prefix}_macs2_peaks.narrowPeak -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
  grep 'mapped (' $sampleFlagstat | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${prefix}", a/\$1}' | cat $fripScoreHeader - > ${prefix}_macs2_peaks.FRiP_mqc.tsv
  """
}


/*
 * GENRICH
 */

process genrich {
  tag "${prefix}"
  label 'genrich'
  label 'medCpu'
  label 'highMem'
  publishDir path: "${params.outDir}/peakCalling/genrich", mode: 'copy',
    saveAs: { filename -> if (filename.endsWith(".tsv")) "stats/$filename"
            else filename }

  when:
  !params.skipPeakCalling && params.caller =~ "genrich"

  input:
  set val(prefix), file(bam), file(sampleFlagstat) from chBamsGenrich.join(chFilteredGenrichFlagstat)
  file peakCountHeader from chPeakCountHeaderGenrich.collect()
  file fripScoreHeader from chFripScoreHeaderGenrich.collect()

  output:
  set val(prefix), file("*Genrich_peaks.narrowPeak"),val("Genrich") into chGenrichPeaks,chGenrichPeaksbb
  file "*_mqc.tsv" into chGenrichMqc
  file("v_genrich.txt") into chGenrichVersion
  file("v_bedtools.txt") into chGenrichBedtoolsVersion

  script:
  opts=params.tn5sites ? "-j -y -d 100" : ""
  """
  echo \$(Genrich --version) &> v_genrich.txt
  echo \$(bedtools --version 2>&1) &> v_bedtools.txt
  samtools sort -n -@ ${task.cpus} -o ${prefix}_nsorted.bam ${bam[0]} 
  Genrich -t ${prefix}_nsorted.bam -o ${prefix}_Genrich_peaks.narrowPeak -D ${opts}
  sed -i '1 i\\track type=narrowPeak name=\"${prefix}_genrich\" description=\"${prefix}_genrich\" nextItemButton=on' ${prefix}_Genrich_peaks.narrowPeak
  cat ${prefix}_Genrich_peaks.narrowPeak | tail -n +2 | wc -l | awk -v OFS='\t' '{ print "${prefix}", \$1 }' | cat $peakCountHeader - > ${prefix}_genrich_peaks.count_mqc.tsv
  READS_IN_PEAKS=\$(intersectBed -a ${bam[0]} -b ${prefix}_Genrich_peaks.narrowPeak -bed -c -f 0.20 | awk -F '\t' '{sum += \$NF} END {print sum}')
  grep 'mapped (' $sampleFlagstat | awk -v a="\$READS_IN_PEAKS" -v OFS='\t' '{print "${prefix}", a/\$1}' | cat $fripScoreHeader - > ${prefix}_genrich_peaks.FRiP_mqc.tsv
  """
}

chMacsPeaks
  .concat(chGenrichPeaks)
  .dump(tag : 'peaks')
  .into{ chPeaksHomer; chPeakQC }

process narrowPeaksToBigBed {
  tag "${prefix}"
  label 'bedtobigbed'
  label 'medCpu'
  label 'medMem'
  publishDir path: "${params.outDir}/peakCalling", mode: 'copy',
    saveAs: { filename -> if (filename.endsWith("Genrich_peaks.narrowPeak.bb")) "genrich/$filename"
            else if (filename.endsWith("macs2_peaks.narrowPeak.bb")) "macs2/$filename"
            else filename }

  input:
  set val(prefix),file(peakfile),val(caller) from chMacsPeaksbb.concat(chGenrichPeaksbb).dump(tag : 'bb')
  file(chrSizes) from chChromSize.collect()

  output:
  set val(prefix), file("*peaks.narrowPeak.bb"),val(caller) into chbigBeds

  script:
  """
  awk -v OFS='\t' '{\$5=\$5>1000?1000:\$5} {print}' ${peakfile} | tail -n +2 | cut -f 1-6 | sort -k 1,1 -k2,2n > ${prefix}_${caller}.tmp
  bedClip -truncate ${prefix}_${caller}.tmp ${chrSizes} ${prefix}_${caller}.bedclip
  bedToBigBed -type=bed6 ${prefix}_${caller}.bedclip ${chrSizes} ${prefix}_${caller}_peaks.narrowPeak.bb
  """
}


/************************************
 * Peaks Annotation
 */

process peakAnnoHomer{
  tag "${sample}"
  label 'homer'
  label 'medCpu'
  label 'medMem'
  publishDir path: "${params.outDir}/peakCalling/annotation/", mode: 'copy'

  when:
  !params.skipPeakAnno

  input:
  set val(sample), file (peakfile), val(caller) from chPeaksHomer
  file gtfFile from chGtfHomer.collect()
  file fastaFile from chFastaHomer.collect()

  output:
  file "*.txt" into chHomerMqc

  script:
  """
  annotatePeaks.pl $peakfile \\
        $fastaFile \\
        -gtf $gtfFile \\
        -cpu ${task.cpus} \\
        > ${sample}_${caller}_annotated_peaks.txt
  """
}


/*
 * Peak calling & annotation QC
 */

process peakQC{
  label 'r'
  label 'medCpu'
  label 'medMem'
  publishDir "${params.outDir}/peakCalling/QC/", mode: 'copy'
  errorStrategy 'ignore'

  when:
  !params.skipPeakAnno

  input:
  file peaks from chPeakQC.collect{ it[-2] }
  file annotations from chHomerMqc.collect().ifEmpty([])
  file peakHeader from chPeakAnnotationHeader.ifEmpty([])

  output:
  file("*.{txt,pdf}") into chMacsQcOutput
  file("*.tsv") optional true into chPeakMqc
  file("v_R.txt") into chRVersion

  script:
  """
  echo \$(R --version 2>&1) &> v_R.txt 
  ${baseDir}/bin/plot_macs_qc.r \\
    -i ${peaks.join(',')} \\
    -s ${peaks.join(',').replaceAll("_peaks.narrowPeak","")} \\
    -o ./ \\
    -p peak
  plot_homer_annotatepeaks.r \\
    -i ${annotations.join(',')} \\
    -s ${annotations.join(',').replaceAll("_annotated_peaks.txt","")} \\
    -o ./ \\
    -p annotatePeaks
  cat $peakHeader annotatePeaks.summary.txt > annotatedPeaks.summary_mqc.tsv
  """
}
    

/*
 * MultiQC
 */

process getSoftwareVersions{
  label 'python'
  label 'minCpu'
  label 'minMem'
  publishDir path: "${params.outDir}/softwareVersions", mode: "copy"

  when:
  !params.skipSoftVersions

  input:
  file 'v_fastqc.txt' from chFastqcVersion.first().ifEmpty([])
  file 'v_bwa.txt' from chBwaVersion.first().ifEmpty([])
  file 'v_bowtie2.txt' from chBowtie2Version.first().ifEmpty([])
  file 'v_samtools.txt' from chSamtoolsVersionBamSort.concat(chSamtoolsVersionBamFiltering).first().ifEmpty([])
  file 'v_picard.txt' from chPicardVersion.first().ifEmpty([])
  file 'v_macs2.txt' from chMacsVersion.first().ifEmpty([])
  file 'v_genrich.txt' from chGenrichVersion.first().ifEmpty([])
  file 'v_preseq.txt' from chPreseqVersion.first().ifEmpty([])
  file 'v_deeptools.txt' from chDeeptoolsVersion.first().ifEmpty([])
  file 'v_bedtools.txt' from chGenrichBedtoolsVersion.concat(chMacsBedtoolsVersion).first().ifEmpty([])
  file 'v_r.txt' from chRVersion.ifEmpty([])

  output:
  file 'software_versions_mqc.yaml' into softwareVersionsYaml

  script:
  """
  echo $workflow.manifest.version &> v_pipeline.txt
  echo $workflow.nextflow.version &> v_nextflow.txt
  scrape_software_versions.py &> software_versions_mqc.yaml
  """
}

process workflowSummaryMqc {
  when:
  !params.skipMultiQC

  output:
  file 'workflow_summary_mqc.yaml' into workflowSummaryYaml

  exec:
  def yaml_file = task.workDir.resolve('workflow_summary_mqc.yaml')
  yaml_file.text  = """
  id: 'summary'
  description: " - this information is collected when the pipeline is started."
  section_name: 'Workflow Summary'
  section_href: 'https://gitlab.curie.fr/data-analysis/atac-seq'
  plot_type: 'html'
  data: |
        <dl class=\"dl-horizontal\">
  ${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
  """.stripIndent()
}

process multiqc {
  label 'multiqc'
  label 'lowCpu'
  label 'lowMem'
  publishDir "${params.outDir}/MultiQC", mode: 'copy'

  when:
  !params.skipMultiQC && !params.inputBam

  input:
  file (splan) from chSplan.collect()
  file metadata from chMetadata.ifEmpty([])
  file multiqcConfig from chMultiqcConfig.ifEmpty([])
  file ('softwareVersions/*') from softwareVersionsYaml.collect().ifEmpty([])
  file ('summary/*') from workflowSummaryYaml.collect()
  file ('fastqc/*') from chFastqcMqc.collect().ifEmpty([])
  file ('mapping/*') from chMappingMqc.collect().ifEmpty([])
  file ('mapping/*') from chMarkedPicstats.collect().ifEmpty([])
  file ('mapping/stats/*') from chRawStats.collect().ifEmpty([])
  file ('mapping/stats/*') from chFilteredStats.collect().ifEmpty([])
  file ('mapping/*') from chStatsMqc.collect().ifEmpty([])
  file ('preseq/*') from chPreseqStats.collect().ifEmpty([])
  file ('fragSize/picard/*') from chFragmentsSize.collect().ifEmpty([])
  file ('deepTools/*') from chDeeptoolsSingleMqc.collect().ifEmpty([])
  file ("deepTools/*") from chDeeptoolsFingerprintMqc.collect().ifEmpty([])
  file ('peakCalling/*') from chMacsOutput.collect().ifEmpty([])
  file ('peakCalling/*') from chMacsMqc.collect().ifEmpty([])
  file ('peakCalling/*') from chGenrichMqc.collect().ifEmpty([])
  file('peakQC/*') from chPeakMqc.collect().ifEmpty([])
  
  output:
  file splan
  file ("*_report.html") into multiqc_report
  file ("*_data")

  script:
  rtitle = customRunName ? "--title \"$customRunName\"" : ''
  rfilename = customRunName ? "--filename " + customRunName + "_atacseq_report" : "--filename atacseq_report" 
  metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
  isPE = params.singleEnd ? "" : "-p"
  designOpts= params.design ? "-d ${params.design}" : ""
  modules_list = "-m custom_content -m fastqc -m bowtie2 -m preseq -m picard -m deeptools -m macs2 -m homer"
  """
  stats2multiqc.sh -s ${splan} ${designOpts} -a ${params.aligner} -m ${params.mitoName} ${isPE}
  medianReadNb="\$(sort -t, -k3,3n mqc.stats | awk -F, '{a[i++]=\$3;} END{x=int((i+1)/2); if (x<(i+1)/2) printf "%.0f", (a[x-1]+a[x])/2; else printf "%.0f",a[x-1];}')"
  mqc_header.py --splan ${splan} --name "ATAC-seq" --version ${workflow.manifest.version} ${metadataOpts} --nbreads \${medianReadNb} > multiqc-config-header.yaml
  multiqc . -f $rtitle $rfilename -c multiqc-config-header.yaml -c $multiqcConfig $modules_list
  """
}

/*
 * Sub-routine
 */
process outputDocumentation {
    label 'python'
    label 'lowCpu'
    label 'lowMem'
    publishDir "${params.summaryDir}/", mode: 'copy'

    input:
    file output_docs from chOutputDocs
    file images from chOutputDocsImages

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

workflow.onComplete {

  /*pipeline_report.html*/

  // pipelineReport.html
  def reportFields = [:]
  reportFields['pipeline'] = workflow.manifest.name
  reportFields['version'] = workflow.manifest.version
  reportFields['runName'] = customRunName ?: workflow.runName
  reportFields['success'] = workflow.success
  reportFields['dateComplete'] = workflow.complete
  reportFields['duration'] = workflow.duration
  reportFields['exitStatus'] = workflow.exitStatus
  reportFields['errorMessage'] = (workflow.errorMessage ?: 'None')
  reportFields['errorReport'] = (workflow.errorReport ?: 'None')
  reportFields['commandLine'] = workflow.commandLine
  reportFields['projectDir'] = workflow.projectDir
  reportFields['summary'] = summary
  reportFields['summary']['Date Started'] = workflow.start
  reportFields['summary']['Date Completed'] = workflow.complete
  reportFields['summary']['Pipeline script file path'] = workflow.scriptFile
  reportFields['summary']['Pipeline script hash ID'] = workflow.scriptId
  if(workflow.repository) reportFields['summary']['Pipeline repository Git URL'] = workflow.repository
  if(workflow.commitId) reportFields['summary']['Pipeline repository Git Commit'] = workflow.commitId
  if(workflow.revision) reportFields['summary']['Pipeline Git branch/tag'] = workflow.revision

  // Render the TXT template
  def engine = new groovy.text.GStringTemplateEngine()
  def tf = new File("$baseDir/assets/onCompleteTemplate.txt")
  def txtTemplate = engine.createTemplate(tf).make(reportFields)
  def reportTxt = txtTemplate.toString()

  // Render the HTML template
  def hf = new File("$baseDir/assets/onCompleteTemplate.html")
  def htmlTemplate = engine.createTemplate(hf).make(reportFields)
  def reportHtml = htmlTemplate.toString()

  // Write summary HTML to a file
  def outputSummaryDir = new File( "${params.summaryDir}/" )
  if( !outputSummaryDir.exists() ) {
    outputSummaryDir.mkdirs()
  }
  def outputHtmlFile = new File( outputSummaryDir, "pipelineReport.html" )
  outputHtmlFile.withWriter { w -> w << reportHtml }
  def outputTxtFile = new File( outputSummaryDir, "pipelineReport.txt" )
  outputTxtFile.withWriter { w -> w << reportTxt }

  // onComplete file
  File woc = new File("${params.outDir}/workflowOnComplete.txt")
  Map endSummary = [:]
  endSummary['Completed on'] = workflow.complete
  endSummary['Duration']     = workflow.duration
  endSummary['Success']      = workflow.success
  endSummary['exit status']  = workflow.exitStatus
  endSummary['Error report'] = workflow.errorReport ?: '-'
  String endWfSummary = endSummary.collect { k,v -> "${k.padRight(30, '.')}: $v" }.join("\n")
  println endWfSummary
  String execInfo = "Execution summary\n${endWfSummary}\n"
  woc.write(execInfo)

  // final logs
  if(workflow.success){
      log.info "Pipeline Complete"
  }else{
    log.info "FAILED: $workflow.runName"
  }
}
