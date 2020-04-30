pipedir=$(cd $(dirname $0); pwd)      # path  of our pipeline
working_dic=.
method=default
index_dic=
thread=8

FastViromeExplorer=`grep 'FastViromeExplorer' ${pipedir}/path.txt|awk -F '=' '{print $2}'`
picard=`grep 'picard' ${pipedir}/path.txt|awk -F '=' '{print $2}'`

while getopts ":1:2:o:t:h:l:k:r:i:m:" opt
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
        l)
		    virus_list="$OPTARG"
		    ;;
		    k)
		    kallisto_list="$OPTARG"
		    ;;
        r)
        file_db="$OPTARG"
        ;;
        i)
        index_dic="$OPTARG"
        ;;
		    m)
        method="$OPTARG"
        ;;
		    ?)
        echo "help"
        echo "requested parameter:"
		    echo "-1       input .fastq file or .fastq.gz file for paired-end reads sequence 1"
        echo "-2       input .fastq file or .fastq.gz file for paired-end reads sequence 2"
		    echo "-h       host genome hisat2 index"
        echo "-l 			 the virus list"
		    echo "-k       kallisto index"
        echo "-r       reference viral genome"
        echo "optional parameter:"
        echo "-o       directory for output files; default is current path"
		    echo "-t       number of threads used; default 1"
        echo "-i       path to store single viral genome reference; default will auto create an file refs in working path"
        echo "-m       assembly method used: default or kallisto or bwa"
		    exit 1;;
	  esac
done
if [ -z ${input2}]
then
    bash ${pipedir}/preprocess.sh -1 ${input1} -o ${working_dic} -t ${thread} -h ${host_ref}
else
    bash ${pipedir}/preprocess.sh -1 ${input1} -2 ${input2} -o ${working_dic} -t ${thread} -h ${host_ref}
fi
label=${input1}
name=(${label//./ })
label=${name##*/}
sample_prefix=${label%%_*}
bash ${pipedir}/virus_identification.sh -p ${sample_prefix} -o ${working_dic} -l ${virus_list} -k ${kallisto_list} -t ${thread}
time bash ${pipedir}/virus_assembly.sh -p ${sample_prefix} -r ${file_db} -o ${working_dic} -t ${thread}
echo 'time default'
echo ${sample_prefix}
bash ${pipedir}/report.sh -p ${sample_prefix} -o ${working_dic} -m ${method}

#time bash ${pipedir}/virus_assembly_kallisto.sh -p ${sample_prefix} -r ${file_db} -o ${working_dic} -t ${thread}
#bash ${pipedir}/report.sh -p ${sample_prefix} -o ${working_dic} -m kallisto
#echo 'time kallisto'
#echo ${sample_prefix}