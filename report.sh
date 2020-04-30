pipedir=$(cd $(dirname $0); pwd)      # path  of our pipeline
working_dic=.
method=default
while getopts ":p:o:m:" opt
do
    case $opt in
            p)
            sample_prefix=$OPTARG
            ;;
            o)
            working_dic="$OPTARG"
            ;;
            m)
            method="$OPTARG"
            ;;
            ?)
            echo "help"
            echo "requested parameters:"
            echo "-p       prefix of input .fq files in preprocess output"
            echo "-m       assembly method used: default or kallisto or bwa;"
            echo "optional parameters:"
            echo "-o       directory for output files; default is current path"
            exit 1;;
      esac
done

##################get sample name from input########################
preprocess=${working_dic}/${sample_prefix}/preprocess
report=${working_dic}/${sample_prefix}/report_${method}
assembly=${working_dic}/${sample_prefix}/assembly_${method}
fve=${working_dic}/${sample_prefix}/fve
################################### copy ###############################

if [ -f ${preprocess}/unmapped_2.fq ];
then
    echo -e "Samples\tReads_Total\tRead_QC_Removed\tReads_mapped" >${report}/reads.qc.summary.txt
    c=`expr 0`
    for count in `ls ${preprocess}/*.sort.bam`; do
        c=`expr $c + 1`
        cp ${preprocess}/${c}.${sample_prefix}.1/*/Images/per_base_sequence_content.png ${report}/${c}-1-per-base-sequence-content.png
        cp ${preprocess}/${c}.${sample_prefix}.1/*/Images/per_sequence_gc_content.png ${report}/${c}-1-per-sequence-gc-content.png
        cp ${preprocess}/${c}.${sample_prefix}.1/*/Images/per_sequence_quality.png ${report}/${c}-1-per-sequence-quality.png
        cp ${preprocess}/${c}.${sample_prefix}.1/*/Images/per_base_quality.png ${report}/${c}-1-per-base-quality.png
        cp ${preprocess}/${c}.${sample_prefix}.2/*/Images/per_base_sequence_content.png ${report}/${c}-2-per-base-sequence-content.png
        cp ${preprocess}/${c}.${sample_prefix}.2/*/Images/per_sequence_gc_content.png ${report}/${c}-2-per-sequence-gc-content.png
        cp ${preprocess}/${c}.${sample_prefix}.2/*/Images/per_sequence_quality.png ${report}/${c}-2-per-sequence-quality.png
        cp ${preprocess}/${c}.${sample_prefix}.2/*/Images/per_base_quality.png ${report}/${c}-2-per-base-quality.png
        
        cp ${preprocess}/${c}.out1.fqcheck ${report}/${c}.out1.fqcheck
        cp ${preprocess}/${c}.out2.fqcheck ${report}/${c}.out2.fqcheck
        cp ${preprocess}/${c}.summary.txt ${report}/${c}.summary.txt
        cp ${fve}/FastViromeExplorer-final-sorted-abundance.tsv ${report}/FastViromeExplorer-final-sorted-abundance.tsv
        cp ${fve}/abundance.tsv ${report}/abundance.tsv
        
        ############################### sort data ###########################
        total_1=`sed -n 's/^\([0-9]*\) sequences processed in total*/\1/p' ${preprocess}/${c}.${sample_prefix}.1.txt`
        removed_1=`sed -n 's/^Number of sequence pairs removed because at least one read was shorter than the length cutoff (20 bp): \([0-9]*\) ([0-9|.|%]*).*/\1/p' ${preprocess}/${c}.${sample_prefix}.1.txt`
        if [ ! $removed_1 ];then
            removed_1=0
        fi
        total_2=`sed -n 's/^\([0-9]*\) sequences processed in total*/\1/p' ${preprocess}/${c}.${sample_prefix}.2.txt`
        removed_2=`sed -n 's/^Number of sequence pairs removed because at least one read was shorter than the length cutoff (20 bp): \([0-9]*\) ([0-9|.|%]*).*/\1/p' ${preprocess}/${c}.${sample_prefix}.2.txt`
        if [ ! $removed_2 ];then
            removed_2=0
        fi
        reads_total=$total_1
        [ $reads_total -gt $total_2 ]&&reads_total=$total_2
        removed=$removed_1
        [ $removed_2 -gt $removed ]&&removed=$removed_2
        reads_trimmed=`expr $reads_total - $removed`
        
        reads_clean=`sed -n 's/^\([0-9]*\) reads;.*/\1/p' ${report}/${c}.summary.txt`
        rate_overall=`sed -n 's/^\([0-9|.]*\)% overall alignment rate/\1/p' ${report}/${c}.summary.txt`
        reads_mapped=$(echo ${reads_clean}*${rate_overall}*0.01|bc)
        
        echo -e "${c}.${sample_prefix}\t${reads_total}\t${removed}\t${reads_mapped}" >>${report}/reads.qc.summary.txt
    done
    python3 ${pipedir}/pipeline.virus.distribution.py ${report} ${sample_prefix} ${c} 1
    
    
    ############################ R report ###############################
    Rscript ${pipedir}/pipeline.plot.R ${report} ${sample_prefix}
    cp ${pipedir}/readsQCAnno.csv ${report}/readsQCAnno.csv
    cp ${pipedir}/covVirAnno.csv ${report}/covVirAnno.csv
    cp ${pipedir}/congtigAnno.csv ${report}/congtigAnno.csv
    cp ${pipedir}/qc.paired.brew ${report}/qc.brew
    cp ${pipedir}/title.brew ${report}/title.brew
    cp ${pipedir}/section2.brew ${report}/section2.brew
    cp ${pipedir}/end.brew ${report}/end.brew
    Rscript ${pipedir}/pipeline.report.R ${report} ${sample_prefix} ${c}
    #rm ${report}/*.brew
else
    echo single
    echo -e "Samples\tReads_Total\tRead_QC_Removed\tReads_mapped" >${report}/reads.qc.summary.txt
    c=`expr 0`
    for count in `ls ${preprocess}/*.sort.bam`; do
        c=`expr $c + 1`
        cp ${preprocess}/${c}.${sample_prefix}/*/Images/per_base_sequence_content.png ${report}/${c}-per-base-sequence-content.png
        cp ${preprocess}/${c}.${sample_prefix}/*/Images/per_sequence_gc_content.png ${report}/${c}-per-sequence-gc-content.png
        cp ${preprocess}/${c}.${sample_prefix}/*/Images/per_sequence_quality.png ${report}/${c}-per-sequence-quality.png
        cp ${preprocess}/${c}.${sample_prefix}/*/Images/per_base_quality.png ${report}/${c}-per-base-quality.png
        
        cp ${preprocess}/${c}.out.fqcheck ${report}/${c}.out.fqcheck
        cp ${preprocess}/${c}.summary.txt ${report}/${c}.summary.txt
        cp ${fve}/FastViromeExplorer-final-sorted-abundance.tsv ${report}/FastViromeExplorer-final-sorted-abundance.tsv
        cp ${fve}/abundance.tsv ${report}/abundance.tsv
        
        ############################### sort data ###########################
        total_1=`sed -n 's/^\([0-9]*\) sequences processed in total*/\1/p' ${preprocess}/${c}.${sample_prefix}.txt`
        removed_1=`sed -n 's/^Number of sequence pairs removed because at least one read was shorter than the length cutoff (20 bp): \([0-9]*\) ([0-9|.|%]*).*/\1/p' ${preprocess}/${c}.${sample_prefix}.txt`
        if [ ! $removed_1 ];then
            removed_1=0
        fi
        
        reads_total=$total_1
        removed=$removed_1
        reads_trimmed=`expr $reads_total - $removed`
        
        reads_clean=`sed -n 's/^\([0-9]*\) reads;.*/\1/p' ${report}/${c}.summary.txt`
        rate_overall=`sed -n 's/^\([0-9|.]*\)% overall alignment rate/\1/p' ${report}/${c}.summary.txt`
        reads_mapped=$(echo ${reads_clean}*${rate_overall}*0.01|bc)
        
        echo -e "${c}.${sample_prefix}\t${reads_total}\t${removed}\t${reads_mapped}" >>${report}/reads.qc.summary.txt
    done
    python3 ${pipedir}/pipeline.virus.distribution.py ${report} ${sample_prefix} ${c} 0
    
    
    ############################ R report ###############################
    Rscript ${pipedir}/pipeline.plot.R ${report} ${sample_prefix}
    cp ${pipedir}/readsQCAnno.csv ${report}/readsQCAnno.csv
    cp ${pipedir}/covVirAnno.csv ${report}/covVirAnno.csv
    cp ${pipedir}/congtigAnno.csv ${report}/congtigAnno.csv
    cp ${pipedir}/qc.single.brew ${report}/qc.brew
    cp ${pipedir}/title.brew ${report}/title.brew
    cp ${pipedir}/section2.brew ${report}/section2.brew
    cp ${pipedir}/end.brew ${report}/end.brew
    Rscript ${pipedir}/pipeline.report.R ${report} ${sample_prefix} ${c}
    #rm ${report}/*.brew

fi