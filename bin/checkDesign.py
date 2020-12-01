#!/usr/bin/env python

import argparse
import csv
import sys
import re
import os

def argsParse():
    """
    Parsing input & outputs CSV files. Also takes in a boolean to indicate if
    the raw reads are single-end or paired-end
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--design", dest="design", help="Design file (csv)", default=None)
    parser.add_argument("-s", "--sampleplan", dest="sampleplan", help="SamplePlan file (csv)")
    parser.add_argument("--singleEnd", help="Specify that input reads are single-end", action="store_true")
    parser.add_argument("--baseDir", help="Base dir if needed", default=".")
    parser.add_argument("--bam", help="Specify that input files are BAM files", action="store_true")
    parser.add_argument("-o","--output", dest ="output" , help="Output design file", default=None)
    args = parser.parse_args()
    inputDesign = args.design
    outputDesign = args.output
    inputData = args.sampleplan
    singleEnd = args.singleEnd
    baseDir = args.baseDir
    inputBam = args.bam
    return inputDesign, outputDesign, inputData, singleEnd, baseDir, inputBam


def check_designs(inputDesign, outputDesign,  inputData, isSingleEnd, baseDir, isInputBam):
    dict_design = {
        'SAMPLEID': [],
        'SAMPLENAME': [],
        'GROUP': [],
        'REPLICATE': []
    }
    if (isSingleEnd):
        dict_reads = {
            'SAMPLEID': [],
            'SAMPLENAME': [],
            'FASTQR1': [],
        }
    else:
        dict_reads = {
            'SAMPLEID': [],
            'SAMPLENAME': [],
            'FASTQR1': [],
            'FASTQR2': []
        }

    ### Checks for design file
    if inputDesign is not None:
        with open(inputDesign, 'r') as designFile:
            lines = csv.reader(designFile)
            header = next(lines)
            for i in range(0, len(header)):
                try:
                    if not header[i] == [*dict_design][i]:
                        raise()
                except:
                    print('Design file columns are not valid, should be : {}'
                          .format([*dict_design]))
                    sys.exit(1)
            # Fill dict to check all input design data
            for sample in lines:
                dict_design['SAMPLEID'].append(sample[0])
                dict_design['SAMPLENAME'].append(sample[1])
                dict_design['GROUP'].append(sample[2])
                dict_design['REPLICATE'].append(sample[3])
            # Check if samples and controls are correctly separated
            #for ID in dict_design['CONTROLID']:
            #    if ID in dict_design['SAMPLEID']:
            #        print('The sample {} is both qualified as a control and a sample'
            #            .format(ID))
            #        sys.exit(1)
            #print("[DESIGN] Check samples/controls IDs ... ok")
            # Check if peaktypes are correct for every sample
            #peaktype_list = ['sharp', 'broad', 'very-broad']
            #index = 0
            #for peaktype in dict_design['PEAKTYPE']:
            #    if not peaktype in peaktype_list:
            #        print('Peaktype for {} is invalid, can be : {}'
            #            .format(dict_design['SAMPLEID'][index], 
            #            ', '.join(peaktype_list)))
            #    index += 1
            #print("[DESIGN] Check peak type information ... ok")

    ### Checks for sampleplan file
    with open(inputData, 'r') as dataFile:
        lines = csv.reader(dataFile)
        # Fill dict to check all input sample data
        for sample in lines:
            dict_reads['SAMPLEID'].append(sample[0])
            dict_reads['SAMPLENAME'].append(sample[1])
            if sample[2][0] != '/':
                readfile = baseDir + '/' + sample[2]
                dict_reads['FASTQR1'].append(readfile)
            else:
                dict_reads['FASTQR1'].append(sample[2])
            if not isSingleEnd:
                if sample[3][0] != '/':
                    readfile = baseDir + '/' + sample[3]
                    dict_reads['FASTQR2'].append(readfile)
                else:
                    dict_reads['FASTQR2'].append(sample[3])

        if inputDesign is not None:
            # Check if there is a missing ID in the design or sample file
            for ID in dict_reads['SAMPLEID']:
                if not (ID in dict_design['SAMPLEID']):
                    print('Sample plan contains an ID that is not in the '
                        'design file ({})'.format(ID))
                    sys.exit(1)
        # Check paths to files
        for samplePath in dict_reads['FASTQR1']:
            if not os.path.exists(samplePath):
                print('File not found {}'
                      .format(os.path.basename(samplePath)))
                sys.exit(1)
        # Check file extensions to match fastq or sam/bam ones
            if(isInputBam == False):
                if not ((samplePath.endswith('fq.gz'))
                    or (samplePath.endswith('fastq.gz'))
                    or (samplePath.endswith('fastq'))
                    or (samplePath.endswith('fq'))):
                    print('File not found {}'
                          .format(os.path.basename(samplePath)))
                    sys.exit(1)
            else:
                if not  ((samplePath.endswith('sam'))
                    or  (samplePath.endswith('bam'))):
                    print('The file {} is not a sam/bam file'
                          .format(os.path.basename(samplePath)))
                    sys.exit(1)
        print("[SAMPLEPLAN] Check file paths ... ok")

    ###### ADDED CHECK UPS ###########
    ### SamplePlan
    ## CHECK IF DATA IS PAIRED-END OR SINGLE-END AND NOT A MIXTURE
    ############## PENSER A ADAPTER CES LIGNES AU SAMPLEPLAN ET NON DESIGN FILE ####################
    #if min(numColList) != max(numColList):
    #    print("{}: Mixture of paired-end and single-end reads!".format(ERROR_STR))
    #    sys.exit(1)


    ### Design
    ## CHECK HEADER
    ERROR_STR = 'ERROR: Please check design file'
    HEADER = ['SAMPLEID', 'SAMPLENAME', 'GROUP', 'REPLICATE']
    fin = open(inputDesign,'r')
    header = fin.readline().strip().split(',')
    if header != HEADER:
        print("{} header: {} != {}".format(ERROR_STR,','.join(header),','.join(HEADER)))
        sys.exit(1)


    numColList = []
    groupRepDict = {}
    while True:
        line = fin.readline()
        if line:
            lspl = [x.strip() for x in line.strip().split(',') if x]

            ## CHECK VALID NUMBER OF COLUMNS PER SAMPLE
            numCols = len(lspl)
            #if numCols not in [3,4]:
            if numCols not in [4]:
                print("{}: Invalid number of columns (Design file must contain 4 comma separated columns, with the following information SAMPLEID,SAMPLENAME,GROUP,REPLICATE)!\nLine: '{}'".format(ERROR_STR,line.strip()))
                sys.exit(1)
            numColList.append(numCols)

            ## CHECK GROUP COLUMN HAS NO SPACES
            #group,replicate,fastQFiles = lspl[0],lspl[1],lspl[2:]
            SAMPLEID,SAMPLENAME,GROUP,REPLICATE = lspl[0],lspl[1],lspl[2],lspl[3]
            if GROUP.find(' ') != -1:
                print("{}: Group id contains spaces!\nLine: '{}'".format(ERROR_STR,line.strip()))
                sys.exit(1)
            ## CHECK SAMPLEID COLUMN HAS NO SPACES
            if SAMPLEID.find(' ') != -1:
                print("{}: sample id contains spaces!\nLine: '{}'".format(ERROR_STR,line.strip()))
                sys.exit(1)
            ## CHECK SAMPLENAME COLUMN HAS NO SPACES
            if SAMPLENAME.find(' ') != -1:
                print("{}: sample names contains spaces!\nLine: '{}'".format(ERROR_STR,line.strip()))
                sys.exit(1)

            ## CHECK REPLICATE COLUMN IS INTEGER
            if not REPLICATE.isdigit():
                print("{}: Replicate id not an integer!\nLine: '{}'".format(ERROR_STR,line.strip()))
                sys.exit(1)

            ## CREATE GROUP MAPPING DICT = {GROUP_ID: {REPLICATE_ID:[[FASTQ_FILES]]}
            REPLICATE = int(REPLICATE)
            if GROUP not in groupRepDict:
                groupRepDict[GROUP] = {}
            #if REPLICATE not in groupRepDict[GROUP]:
            #    groupRepDict[GROUP][REPLICATE] = []
            #groupRepDict[GROUP][REPLICATE].append(fastQFiles)
            if REPLICATE not in groupRepDict[GROUP]:
                groupRepDict[GROUP][REPLICATE] = []
            groupRepDict[GROUP][REPLICATE].append(SAMPLENAME)


        else:
            fin.close()
            break

    ## CHECK IF MULTIPLE GROUPS EXIST
    multiGroups = False
    if len(groupRepDict) > 1:
        multiGroups = True

    ## WRITE TO FILE
    numRepList = []
    fout = open(outputDesign,'w')
    #fout.write(','.join(['sample_id','fastq_1','fastq_2']) + '\n')
    fout.write(',' + '\n')
    for GROUP in sorted(groupRepDict.keys()):

        ## CHECK THAT REPLICATE IDS ARE IN FORMAT 1..<NUM_REPLICATES>
        uniq_rep_ids = set(groupRepDict[GROUP].keys())
        if len(uniq_rep_ids) != max(uniq_rep_ids):
            print("{}: Replicate IDs must start with 1..<num_replicates>\nGroup: {}, Replicate IDs: {}".format(ERROR_STR,GROUP,list(uniq_rep_ids)))
            sys.exit(1)
        numRepList.append(max(uniq_rep_ids))

        for REPLICATE in sorted(groupRepDict[GROUP].keys()):
            for idx in range(len(groupRepDict[GROUP][REPLICATE])):
                #fastQFiles = groupRepDict[GROUP][REPLICATE][idx]
                SAMPLENAME = groupRepDict[GROUP][REPLICATE][idx]
                SAMPLEID = "{}_R{}_T{}".format(GROUP,REPLICATE,idx+1)
                #if len(fastQFiles) == 1:
                    #fout.write(','.join([sample_id] + fastQFiles) + ',\n')
                #    fout.write(','.join([SAMPLEID] + SAMPLENAME) + ',\n')
                #else:
                    #fout.write(','.join([sample_id] + fastQFiles) + '\n')
                #    fout.write(','.join([SAMPLEID] + SAMPLENAME) + '\n')
                fout.write(','.join([SAMPLEID] + [SAMPLENAME]) + '\n')
    fout.close()

    ## CHECK IF REPLICATES IN DESIGN
    repsExist = False
    if max(numRepList) != 1:
        repsExist = True

    ## CHECK FOR BALANCED DESIGN ACROSS MULTIPLE GROUPS.
    balancedDesign = False
    if len(set(numRepList)) == 1 and multiGroups and repsExist:
        balancedDesign = True


if __name__ == '__main__':
    inputDesign, outputDesign, inputData, singleEnd, baseDir, inputBam = argsParse()
    check_designs(inputDesign, outputDesign, inputData, singleEnd, baseDir, inputBam)
