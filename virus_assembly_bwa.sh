pipedir=$(cd $(dirname $0); pwd)      # path  of our pipeline
picard=`grep 'picard' ${pipedir}/path.txt|awk -F '=' '{print $2}'`
thread=1
working_dic=.
index_dic=
while getopts ":p:o:r:i:t:" opt
do
    case $opt in
            p)
            sample_prefix=$OPTARG
            ;;
            o)
            working_dic="$OPTARG"
            ;;
            t)
            thread="$OPTARG"
            ;;
            r)
            file_db="$OPTARG"
            ;;
            i)
            index_dic="$OPTARG"
            ;;
            ?)
            echo "help"
            echo "requested parameters:"
            echo "-p       prefix of input .fq files in preprocess output"
            echo "-r       reference viral genome"
            echo "-i       path to store single viral genome reference; default will auto create an file refs in working path"
            echo "optional parameters:"
            echo "-o       directory for output files; default is current path"
            echo "-t       number of threads used; default 1"
            exit 1;;
      esac
done

# make dictory for assembly
report=${working_dic}/${sample_prefix}/report_bwa/          # dictory for report generation outputs
mkdir $report
preprocess=${working_dic}/${sample_prefix}/preprocess       # dictory for preprocess step outputs
fve=${working_dic}/${sample_prefix}/fve                     # dictory for virus_identification step outputs
assembly=${working_dic}/${sample_prefix}/assembly_bwa       # dictory for virus assembly step outputs
mkdir ${assembly}
refs=${working_dic}/${sample_prefix}/${index_dic}/refs     # create dictory if not exist to store each identified virus from virus identification step
mkdir ${refs}

################## virus assembly ###################################
Accs=`awk '{print $1}' ${fve}/FastViromeExplorer-final-sorted-abundance.tsv|grep -v '#'`
for acc in ${Accs[@]};do
    ref=${refs}/${acc}.fa
    ###add index for new virus####
		if [ ! -f ${ref} ]
    then
		    echo "$acc" >${assembly}/temp.${acc}.list.txt
      	perl -e 'open LST, shift;while(<LST>){chomp;$hash{$_}=1} close LST;$/=">";<>;while(<>){chomp;$id=(split)[0];if(exists($hash{$id})){print ">",$_;}}' ${assembly}/temp.${acc}.list.txt ${file_db} >${ref}
        bwa index -p ${refs}/${acc} -a is ${ref}
    fi
    
    ### bwa assembly of virome sequence ###
    if [ ! -f ${preprocess}/unmapped_2.fq ];
    then
        echo single
        bwa aln -t ${thread} ${refs}/${acc} ${preprocess}/unmapped.fq > ${assembly}/${acc}.sai
        bwa samse ${refs}/${acc} ${assembly}/${acc}.sai ${preprocess}/unmapped.fq > ${assembly}/${acc}.sam
    else
        echo paired
        bwa aln -t ${thread} ${refs}/${acc} ${preprocess}/unmapped_1.fq > ${assembly}/${acc}_bwa_1.sai
        bwa aln -t ${thread} ${refs}/${acc} ${preprocess}/unmapped_2.fq > ${assembly}/${acc}_bwa_2.sai
        bwa sampe ${refs}/${acc} ${assembly}/${acc}_bwa_1.sai ${assembly}/${acc}_bwa_2.sai ${preprocess}/unmapped_1.fq ${preprocess}/unmapped_2.fq > ${assembly}/${acc}.sam
    fi
    #rm ${assembly}/${acc}_bwa*
    
    
    #remove duplicates and sort
    samtools view -@ ${thread} -F 4 -Sbh ${assembly}/${acc}.sam|samtools sort -@ ${thread} -o ${assembly}/${acc}.bam 
    java -jar ${picard} MarkDuplicates M=${assembly}/${acc}_dupstats REMOVE_DUPLICATES=TRUE I=${assembly}/${acc}.bam O=${assembly}/${acc}_nodup.bam
    samtools index ${assembly}/${acc}_nodup.bam
    #rm ${assembly}/${acc}.sam
    #rm ${assembly}/${acc}.bam

    ### call snvs in virus ####
    rm ${assembly}/${acc}.cov.txt
    rm ${assembly}/${acc}.basecov.txt
    pileup.sh in=${assembly}/${acc}_nodup.bam out=${assembly}/${acc}.cov.txt basecov=${assembly}/${acc}.basecov.txt
    grep -v '#' ${assembly}/${acc}.basecov.txt|awk '{print $1"\t"$2+1"\t"$2+1"\t"$3}' |awk '$4<=3 {print}' >${assembly}/${acc}_basecov_bed 
    bcftools mpileup -f ${refs}/${acc}.fa ${assembly}/${acc}_nodup.bam | bcftools call -mv -Oz -o ${assembly}/${acc}.vcf.gz  
    tabix ${assembly}/${acc}.vcf.gz
    ###produce a fa file of sequence####
    cat ${refs}/${acc}.fa | bcftools consensus --iupac-codes -H 2pIu -m ${assembly}/${acc}_basecov_bed ${assembly}/${acc}.vcf.gz > ${assembly}/${acc}.fa
    
    cp ${assembly}/${acc}.basecov.txt ${report}/${acc}.basecov.txt    
    
done
cat ${assembly}/*.fa >${report}/sequence.fa
cat ${assembly}/*.cov.txt|grep -v '#ID' >${report}/temp.cov.clean.txt
sed -n '1p' ${assembly}/${acc}.cov.txt >${report}/temp.header
cat ${report}/temp.header ${report}/temp.cov.clean.txt >${report}/cov.summary.txt
rm ${report}/temp.cov.clean.txt
rm ${report}/temp.header