pipedic=$(cd $(dirname $0); pwd)      # path  of our pipeline
thread=1
FastViromeExplorer=`grep 'FastViromeExplorer' ${pipedic}/path.txt|awk -F '=' '{print $2}'`
working_dic=.
while getopts ":p:o:l:t:k:" opt
do
    case $opt in
		p)
		sample_prefix=$OPTARG
		;;
		o)
		working_dic=$OPTARG
		;;
		l)
		virus_list="$OPTARG"
		;;
		t)
		thread="$OPTARG"
		;;
		k)
		kallisto_list="$OPTARG"
		;;
    ?)
		echo "help"
    echo "requested parameters:"
		echo "-p 	      	prefix of input .fq files in preprocess output"
		echo "-l 			    the virus list"
		echo "-k        	kallisto index"
		echo "optional parameters:"
    echo "-o       		directory for output files; default is current path"
		echo "-t 			    number of threads used; default 1"
		exit 1;;
	esac
done

preprocess=${working_dic}/${sample_prefix}/preprocess				# path of results from preprocess step
fve=${working_dic}/${sample_prefix}/fve/							      # path of virus identification step, would created automaticly below

################## virus detection #################################
if [ -f ${preprocess}/unmapped_2.fq ];
then
    echo paired
    java ${FastViromeExplorer} -1 ${preprocess}/unmapped_1.fq -2 ${preprocess}/unmapped_2.fq -l ${virus_list} -i ${kallisto_list} -o $fve
else
    java ${FastViromeExplorer} -1 ${preprocess}/unmapped.fq -l ${virus_list} -i ${kallisto_list} -o $fve
fi