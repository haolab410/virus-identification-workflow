working_dir=.
thread=1
while getopts ":1:2:o:t:h:" opt
do
    case $opt in
		    1)
		    input1=$OPTARG
		    ;;
		    2)
		    input2=$OPTARG
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
		    ?)
		    echo "help"
        echo "requested parameter:"
		    echo "-1       input .fastq file or .fastq.gz file for paired-end reads sequence 1"
        echo "-2       input .fastq file or .fastq.gz file for paired-end reads sequence 2"
		    echo "-h       host genome hisat2 index"
        echo "optional parameter:"
        echo "-o       directory for output files; default is current path"
		    echo "-t       number of threads used; default 1"
		    exit 1;;
	  esac
done

# get the sample prefix for naming
label=${input1}
name=(${label//./ })
label=${name##*/}
sample_prefix=${label%%_*}

# make dictory for preprocess
mkdir ${working_dic}/${sample_prefix}			# create new working directory for the sample
preprocess=${working_dic}/${sample_prefix}/preprocess
mkdir ${preprocess}								# create a preprocess file in the working directory for all the results in preprocess step

################## quality control #################################
fastqc -o ${preprocess}/ -t ${thread} ${input1} ${input2}
cd ${preprocess}/
unzip ${sample_prefix}*1*.zip
unzip ${sample_prefix}*2*.zip

FxTools_Linux Fqtools fqcheck -i ${input1} ${input2} -o ${preprocess}/out1 ${preprocess}/out2
trim_galore --fastqc --paired ${input1} ${input2} -o ${preprocess}/ --gzip

################## remove source sequence #################################
hisat2 -x ${host_ref} -q -1 ${preprocess}/${sample_prefix}*1_val_1*.gz -2 ${preprocess}/${sample_prefix}*2_val_2*.gz -p ${thread} -S ${preprocess}/mapped.sam 2> ${preprocess}/summary.txt

samtools view -u -@ ${thread} -f 4 -F 8 ${preprocess}/mapped.sam|samtools sort -@ ${thread} -o ${preprocess}/unmapped_1.sort.bam 					#mapped with R1 but not R2
samtools view -u -@ ${thread} -f 8 -F 4 ${preprocess}/mapped.sam|samtools sort -@ ${thread} -o ${preprocess}/unmapped_2.sort.bam 					#mapped with R2 but not R2
samtools view -u -@ ${thread} -f 12 -F 256 ${preprocess}/mapped.sam|samtools sort -@ ${thread} -o ${preprocess}/unmapped_both.sort.bam 				#unmapped with both R1 and R2
samtools merge ${preprocess}/unmapped.bam ${preprocess}/unmapped_1.sort.bam ${preprocess}/unmapped_2.sort.bam ${preprocess}/unmapped_both.sort.bam  #merge the above three files to unmapped.bam
rm ${preprocess}/*.sort.bam
samtools sort -n ${preprocess}/unmapped.bam -o ${preprocess}/unmapped.sort.bam
rm ${preprocess}/unmapped.bam

bam2fastq ${preprocess}/unmapped.sort.bam -o ${preprocess}/unmapped#.fq
