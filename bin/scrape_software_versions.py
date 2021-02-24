#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re
from os import path

# TODO nf-core: Add additional regexes for new tools in process get_software_versions
regexes = {
    'Pipeline': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
    'BWA': ['v_bwa.txt', r"Version: (\S+)"],
    'Bowtie2': ['v_bowtie2.txt', r"version (\S+)"],
    'samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'bedtools': ['v_bedtools.txt', r"bedtools v(\S+)"],
    'picard': ['v_picard.txt', r"([\d\.]+)-SNAPSHOT"],
    'preseq': ['v_preseq.txt', r"Version: (\S+)"],
    'deeptools': ['v_deeptools.txt', r"bamCoverage (\S+)"],
    'R': ['v_R.txt', r"R version (\S+)"],
    'MACS2': ['v_macs2.txt', r"macs2 (\S+)"],
    'Genrich' : ["v_genrich.txt",r"genrich (\$+)"],
}


results = OrderedDict()
results['Pipeline'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['BWA'] = '<span style="color:#999999;\">N/A</span>'
results['Bowtie2'] = '<span style="color:#999999;\">N/A</span>'
results['samtools'] = '<span style="color:#999999;\">N/A</span>'
results['bedtools'] = '<span style="color:#999999;\">N/A</span>'
results['picard'] = '<span style="color:#999999;\">N/A</span>'
results['preseq'] = '<span style="color:#999999;\">N/A</span>'
results['deeptools'] = '<span style="color:#999999;\">N/A</span>'
results['R'] = '<span style="color:#999999;\">N/A</span>'
results['MACS2'] = '<span style="color:#999999;\">N/A</span>'
results['Genrich'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    if path.isfile(v[0]):
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))

# Dump to YAML
print ('''
id: 'software versions'
section_name: 'Software Versions'
section_href: 'https://gitlab.curie.fr/data-analysis/atac-seq/'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd>{}</dd>".format(k,v))
print ("    </dl>")
