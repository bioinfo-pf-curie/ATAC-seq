version-1.0.4

NEW FEATURES

  - Start analysis from aligned bam files
  - Trim-galore step was added

BUG FIXES

  - Fix bug in clipBed for some genomes
  - Prevent macs2 to export bam files into bed format in the case of the use of the --tn5sites option

***********************************
version-1.0.3

NEW FEATURES

  - add Bombyx Mori V4 annotation

**********************************
version-1.0.2

BUG FIXES

  - Fix bug when no NBR fragments are found

NEW FEATURES

  - add dmelr6.28 genome is now available
  - add plotHeatmap function

***********************************
version-1.0.1

NEW FEATURES

  - new --ignoreBlacklist option to avoid filtering on blacklisted region
  - add --bwaOpts and --bowtie2Opts parameters to change the default mapping options

SIGNIFICANT USER-VISIBLE CHANGES

  - Fix environment issue with tbb and bowtie2
  - Fix a limit in preseq graphical representation in MultiQC and add the median number of reads

***********************************
version-1.0.0

NEW FEATURES

  - First stable of the atac-seq pipeline


