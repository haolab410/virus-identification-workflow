working_dir=.
thread=1
paired='p'

while getopts ":1:2:o:t:h:" opt
do
    case $opt in
		    1)
		    i1=$OPTARG
        ;;
		    o)
    		working_dic="$OPTARG"
		    ;;
		    t)
		    thread="$OPTARG"
		    ;;
        h)
		    host_ref="$OPTARG"
		    ;;
        2)
		    i2=$OPTARG
        ;;
		    ?)
		    echo "help"
        echo "requested parameter:"
		    echo "-1       input .fastq file or .fastq.gz file for reads sequence 1"
        echo "-2       input .fastq file or .fastq.gz file for paired-end reads sequence 2; ignore this if single-end"
		    echo "-h       host genome hisat2 index"
        echo "optional parameter:"
        echo "-o       directory for output files; default is current path"
		    echo "-t       number of threads used; default 1"
		    exit 1;;
	  esac
done

echo ${i2}
if [ -z ${i2}]
then
    echo 'single-end'
    # get the sample prefix for naming
    label=${i1//,/ }
    count=`expr 0`
    for input1 in ${label}; do
        count=`expr $count + 1`
        name=(${input1//.fastq/ })
        name=(${name//.fq/ })
        label=${name##*/}
        sample_prefix=${label%%_*}
        
        # make dictory for preprocess
        mkdir ${working_dic}/${sample_prefix}			# create new working directory for the sample
        preprocess=${working_dic}/${sample_prefix}/preprocess
        mkdir ${preprocess}								# create a preprocess file in the working directory for all the results in preprocess step
        
        fastqc -o ${preprocess}/ -t ${thread} ${input1}
        cd ${preprocess}/
        unzip ${label}*fastqc.zip -d ${count}.${sample_prefix}
        

        FxTools_Linux Fqtools fqcheck -i ${input1} -o ${preprocess}/${count}.out
        trim_galore --fastqc ${input1} -o ${preprocess}/ --gzip --cores 7
        mv ${preprocess}/${label}.*.gz_trimming_report.txt ${preprocess}/${count}.${sample_prefix}.txt
        
        ################## remove source sequence #################################
        hisat2 -x ${host_ref} -q -U ${preprocess}/${label}*trimmed*.gz -p ${thread} -S ${preprocess}/mapped.sam 2> ${preprocess}/${count}.summary.txt
        samtools view -u -@ ${thread} -f 4 ${preprocess}/mapped.sam|samtools sort -@ ${thread} -o ${preprocess}/${count}.sort.bam 		#unmapped reads
        rm ${preprocess}/*.sam
    done
    flag=${preprocess}/2.sort.bam
    if [ ! -f ${flag} ];
    then
        bam2fastq ${preprocess}/1.sort.bam -o ${preprocess}/unmapped.fq     

    else
        samtools merge ${preprocess}/unmapped.bam ${preprocess}/*.sort.bam
        bam2fastq ${preprocess}/unmapped.bam -o ${preprocess}/unmapped#.fq     
    fi

else
    echo 'paired-end'
    label=${i1//,/ }
    count=`expr 0`
    for input1 in ${label}; do
    
        # get the prefix of reads 1
        name=(${input1//.fastq/ })
        name=(${name//.fq/ })
        label1=${name##*/}
        sample_prefix=${label1%%_*}
        # get the prefix of read2
        count=`expr $count + 1`
        input2=`echo $i2|awk -F ',' '{print $"'$count'"}'`
        name=(${input2//.fastq/ })
        name=(${name//.fq/ })
        label2=${name##*/}
        
        
        # make dictory for preprocess
        mkdir ${working_dic}/${sample_prefix}			# create new working directory for the sample
        preprocess=${working_dic}/${sample_prefix}/preprocess
        mkdir ${preprocess}								# create a preprocess file in the working directory for all the results in preprocess step
        
        fastqc -o ${preprocess}/ -t ${thread} ${input1} ${input2}
        cd ${preprocess}/
        unzip ${label1}*fastqc.zip -d ${count}.${sample_prefix}.1
        unzip ${label2}*fastqc.zip -d ${count}.${sample_prefix}.2

        FxTools_Linux Fqtools fqcheck -i ${input1} ${input2} -o ${preprocess}/${count}.out1 ${preprocess}/${count}.out2
        trim_galore --fastqc --paired ${input1} ${input2} -o ${preprocess}/ --gzip
        mv ${preprocess}/${label1}.fq.gz_trimming_report.txt ${preprocess}/${count}.${sample_prefix}.1.txt
        mv ${preprocess}/${label2}.fq.gz_trimming_report.txt ${preprocess}/${count}.${sample_prefix}.2.txt
        
        ################## remove source sequence #################################
        hisat2 -x ${host_ref} -q -1 ${preprocess}/${label1}_val_1*.gz -2 ${preprocess}/${label2}_val_2*.gz -p ${thread} -S ${preprocess}/mapped.sam 2> ${preprocess}/${count}.summary.txt

        samtools view -u -@ ${thread} -f 4 -F 264 ${preprocess}/mapped.sam|samtools sort -@ ${thread} -o ${preprocess}/unmapped_1.sort.bam 					#mapped with R1 but not R2
        samtools view -u -@ ${thread} -f 8 -F 260 ${preprocess}/mapped.sam|samtools sort -@ ${thread} -o ${preprocess}/unmapped_2.sort.bam 					#mapped with R2 but not R2
        samtools view -u -@ ${thread} -f 12 -F 256 ${preprocess}/mapped.sam|samtools sort -@ ${thread} -o ${preprocess}/unmapped_both.sort.bam 				#unmapped with both R1 and R2
        samtools merge ${preprocess}/unmapped.bam ${preprocess}/unmapped_1.sort.bam ${preprocess}/unmapped_2.sort.bam ${preprocess}/unmapped_both.sort.bam  #merge the above three files to unmapped.bam
        rm ${preprocess}/*.sort.bam
        samtools sort -n ${preprocess}/unmapped.bam -o ${preprocess}/${count}.sort.bam
        rm ${preprocess}/unmapped.bam
    done
    flag=${preprocess}/2.sort.bam
    if [ ! -f ${flag} ];
    then
        # only one file
        bam2fastq ${preprocess}/1.sort.bam -o ${preprocess}/unmapped#.fq     
    else
        # multiple file need to be merged
        samtools merge ${preprocess}/unmapped.bam ${preprocess}/*.sort.bam
        bam2fastq ${preprocess}/unmapped.bam -o ${preprocess}/unmapped#.fq    
    fi
    
fi
