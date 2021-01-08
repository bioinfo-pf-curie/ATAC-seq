#!/bin/bash

function usage {
    echo -e "usage : stats2multiqc.sh -s SAMPLE_PLAN -d DESIGN -a ALIGNER [-p][-h]"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo 
    echo "stat2multiqc.sh"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -s SAMPLE_PLAN"
    echo "   -d DESIGN"
    echo "   -a ALIGNER"
    echo "   [-p]: paired-end mode"
    echo "   [-m]: Mitochondrial chromosome name"
    echo "   [-h]: help"
    exit;
}

is_pe=0
mito_name="chrM"
while getopts "s:d:a:m:ph" OPT
do
    case $OPT in
        s) splan=$OPTARG;;
	d) design=$OPTARG;;
	a) aligner=$OPTARG;;
	m) mito_name=$OPTARG;;
	p) is_pe=1;;
	h) help ;;
	\?)
	    echo "Invalid option: -$OPTARG" >&2
	    usage
	    exit 1
	    ;;
	:)
	    echo "Option -$OPTARG requires an argument." >&2
	    usage
	    exit 1
	    ;;
    esac
done

if  [[ -z $splan ]]; then
    usage
    exit
fi

all_samples=$(awk -F, '{print $1}' $splan)

echo -e "Sample_ID,Sample_name,Number_of_reads,Number_of_aligned_reads,Percent_of_aligned_reads,Number_of_mito,Percent_of_mito,Number_of_hq_mapped_reads,Percent_of_hq_mapped_reads,Number_of_lq_mapped_reads,Percent_of_lq_mapped_reads,Number_of_duplicates,Percent_of_duplicates,Number_of_usable_reads,Percent_of_usable_reads,TSS_enrichment,Fraction_of_reads_in_peaks" > mqc.stats

for sample in $all_samples
do
    #SAMPLE NAME
    sname=$(grep "$sample," $splan | awk -F, '{print $2}')

    #ALIGNMENT
    if [ $aligner == "bowtie2" ]; then
	nb_frag=$(grep "reads;" mapping/${sample}_bowtie2.log | sed 's/ .*//')
	if [[ $is_pe == 1 ]]; then
            nb_reads=$(( $nb_frag * 2 ))
	else
            nb_reads=$nb_frag
	fi
    elif [ $aligner == "bwa-mem" ]; then
	# bwa.log file is in reads number (not pairs)
	nb_reads=$(grep 'Total' mapping/${sample}_bwa.log | awk -F "\t" '{print $2}')
	if [[ $is_pe == 1 ]]; then
	    nb_frag=$(( $nb_reads / 2 ))
	else
	    nb_frag=$nb_reads
	fi
	tail -n +3 mapping/${sample}_bwa.log > mapping/${sample}_bwa.mqc
    elif [ $aligner == "star" ]; then
	nb_frag=$(grep "Number of input reads" mapping/${sample}Log.final.out | cut -d"|" -f 2 | sed -e 's/\t//g')
	if [[ $is_pe == 1 ]]; then
            nb_reads=$(( $nb_frag * 2 ))
	else
            nb_reads=$nb_frag
	fi
    fi

    #Mapping stats (always in reads - so must be converted for PE)
    #These statistics are calculated after spike cleaning but before filtering
    nb_mapped=$(awk -F, '$1=="Mapped"{print $2}' mapping/${sample}_mappingstats.mqc)
    nb_mapped_hq=$(awk -F, '$1=="HighQual"{print $2}' mapping/${sample}_mappingstats.mqc)
    nb_mapped_lq=$(awk -F, '$1=="LowQual"{print $2}' mapping/${sample}_mappingstats.mqc)
    perc_mapped=$(echo "${nb_mapped} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    perc_mapped_hq=$(echo "${nb_mapped_hq} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    perc_mapped_lq=$(echo "${nb_mapped_lq} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')

    ##MITO
    if [[ -e mapping/stats/${sample}_raw.idxstats ]]; then
	nb_mito=$(awk -v mt=${mito_name} '$1==mt{print $3}' mapping/stats/${sample}_raw.idxstats)
	perc_mito=$(echo "${nb_mito} ${nb_mapped}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    else
	nb_mito='NA'
	perc_mito='NA'
    fi

    #PICARD
    if [[ -e mapping/${sample}.MarkDuplicates.metrics.txt ]]; then
	nb_dups_pair=$(grep -a2 "## METRICS" mapping/${sample}.MarkDuplicates.metrics.txt | tail -1 | awk -F"\t" '{print $7}')
	nb_dups_single=$(grep -a2 "## METRICS" mapping/${sample}.MarkDuplicates.metrics.txt | tail -1 | awk -F"\t" '{print $6}')
	nb_dups_optical=$(grep -a2 "## METRICS" mapping/${sample}.MarkDuplicates.metrics.txt | tail -1 | awk -F"\t" '{print $8}')
	nb_dups=$(( $nb_dups_pair * 2 + $nb_dups_single + $nb_dups_optical ))
	perc_dups=$(echo "${nb_dups} ${nb_mapped}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    else
	nb_dups='NA'
	perc_dups='NA'
    fi

    #Filtered bam
    if [[ -e mapping/stats/${sample}_filtered.stats ]]; then
	nb_filt=$(grep ^SN mapping/stats/${sample}_filtered.stats | cut -f 2- | grep "reads mapped:" | cut -f 2)
	perc_filt=$(echo "${nb_filt} ${nb_mapped}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    fi

    #TSS Enrichment
    if [[ -e deeptools/${sample}.plotProfile_corrected.tab ]]; then
	tsse=$(awk -F"\t" 'NR==3{mn=mx=$3;for(i=3;i<=NF;i++){if($i>mx){mx=$i}};printf("%.*f",2, mx-mn)}' deeptools/${sample}.plotProfile_corrected.tab)
    else
	tsse='NA'
    fi

    #PeakCalling 
    if [[ ! -z $design ]];
    then
	peaktype=$(grep "^$sample" ${design} | awk -F, '{print $5}') 
	if [ -e peakCalling/${peaktype}/${sample}_peaks.FRiP_mqc.tsv ]; then
	    frip=$(grep "$sample" peakCalling/${peaktype}/${sample}_peaks.FRiP_mqc.tsv | awk '{print $2}')
	else
	    frip='NA'
	fi
    else
	frip='NA'
    fi

    #To file
    echo -e ${sample},${sname},${nb_frag},${nb_mapped},${perc_mapped},${nb_mito},${perc_mito},${nb_mapped_hq},${perc_mapped_hq},${nb_mapped_lq},${perc_mapped_lq},${nb_dups},${perc_dups},${nb_filt},${perc_filt},${tsse},${frip} >> mqc.stats
done

