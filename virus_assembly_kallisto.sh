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
            echo "optional parameters:"
            echo "-i       path to store single viral genome reference; default will auto create an file refs in working path"
            echo "-o       directory for output files; default is current path"
            echo "-t       number of threads used; default 1"
            exit 1;;
      esac
done

# make dictory for assembly
report=${working_dic}/${sample_prefix}/report_kallisto/     # dictory for report generation outputs
mkdir $report
fve=${working_dic}/${sample_prefix}/fve                     # dictory for virus_identification step outputs
assembly=${working_dic}/${sample_prefix}/assembly_kallisto  # dictory for virus assembly step outputs
mkdir ${assembly}
refs=${working_dic}/${sample_prefix}/${index_dic}/refs      # create dictory if not exist to store each identified virus from virus identification step
mkdir ${refs}

# get bam file for virus genome
file_sam=${fve}/FastViromeExplorer-reads-mapped-sorted.sam
samtools view -@ ${thread} -b $file_sam -o ${assembly}/kalli.bam
samtools index -@ ${thread} ${assembly}/kalli.bam
echo 'here'
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
    
    ### extract virome sequence for particular virus id ###
    samtools view -@ ${thread} ${assembly}/kalli.bam ${acc} -o ${assembly}/${acc}_kalli.bam
    # change bam header
    samtools view -H ${assembly}/${acc}_kalli.bam > ${assembly}/temp.sam
    grep 'SO' ${assembly}/temp.sam >${assembly}/temp.ext.header.hd
    grep '@SQ' ${assembly}/temp.sam  |grep $acc >${assembly}/temp.ext.header.acc
    grep '@PG' ${assembly}/temp.sam  >${assembly}/temp.ext.header.pg
    samtools view ${assembly}/${acc}_kalli.bam|grep -v '^@' >${assembly}/temp.ext.reads
    cat ${assembly}/temp.ext.header.hd ${assembly}/temp.ext.header.acc ${assembly}/temp.ext.header.pg ${assembly}/temp.ext.reads >${assembly}/${acc}_kalli.sam
    #rm ${assembly}/temp*
    samtools view -@ ${thread} ${assembly}/${acc}_kalli.sam -o ${assembly}/${acc}_kalli.bam
    #rm ${assembly}/${acc}_kalli.sam
    samtools sort -@ ${thread} ${assembly}/${acc}_kalli.bam -o ${assembly}/${acc}_kalli_sorted.bam
    #rm ${assembly}/${acc}_kalli.bam
    java -jar ${picard} MarkDuplicates M=${assembly}/${acc}_dupstats REMOVE_DUPLICATES=TRUE I=${assembly}/${acc}_kalli_sorted.bam O=${assembly}/${acc}_nodup.bam
    #rm ${assembly}/${acc}_kalli*
    samtools index ${assembly}/${acc}_nodup.bam
    
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