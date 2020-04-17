# virus-identification-workflow

Our virus identification workflow is able to detect viruses from the Next Generation Sequencing raw data. Users can run with raw RNA-seq data towards a final report with data quality, viruses infected and their sequence information. 

# required software:
fastqc; fqcheck 2.05 (https://github.com/BGI-shenzhen/FxTools); 

trim_galore;

hisat2; 

samtools 1.6 and above; 

bam2fastq; 

picard; 

FastViromeExplorer; 

bbmap pileup.sh; 

bwa; 

bcftools 1.8 and above; 

python 3; 

R 3.5 with library brew, ggplot2, geneplotter, optparse.  

You will need to install these softwares and put them in the environment you run.

FastViromeExplorer; picard should be install and users need to put the path of these two software in path.txt.

NCBItax2lin -- can be downloaded from (https://github.com/zyxue/ncbitax2lin), it requires python 2 to run this software. It is used to generate a taxonomy list for later used, and is not included in the workflow. 

# generate required reference:
First, follow the instruction from NCBItax2lin to download taxonomy list lineages-XXXX-XX-XX.csv, this needs to be done in python 2 environment.

    make
    
After make, you will get the file lineages-XXXX-XX-XX.csv

Then, run taxonomy.retriving.sh 

    bash taxonomy.retriving.sh <lineages-XXXX-XX-XX.csv>

to obtain the following files:

  1. virus genome reference 
  2. virus list: viruses-list.txt
  3. kallisto index: ncbi-virus-kallisto-index-k31.idx
# Run pipeline with test data
The test data is from SARS-CoV-2 patients. You can run our pipeline by running pipeline.sh.
        
        bash pipeline.sh -1 test/test_R1.fastq.gz -2 test/test_R2.fastq.gz -o test/ -t 8 -h your_path/grch38/genome -l your_path/viruses-list.txt -k your_path/ncbi-virus-kallisto-index-k31.idx -r your_path/viral.genomic.fna

You can also run step by steps:

   for quality control and remove host sequences:
   
        bash perprocess.sh -1 test/test_R1.fastq.gz -2 test/test_R2.fastq.gz -o test/ -t 8 -h your_path/grch38/genome
   for virus identification:
   
        bash virus_identification.sh -p test -o test/ -l your_path/viruses-list.txt -k your_path/ncbi-virus-kallisto-index-k31.idx -t 8
   for virus assembly with our default method (kallisto + bwa):
   
        bash virus_assembly.sh -p test -r your_path/viral.genomic.fna -o test/ -t 8
   or with kallisto alignment only:
      
        bash virus_assembly_kallisto.sh -p test -r your_path/viral.genomic.fna -o test/ -t 8
   or with bwa alignment only:
   
        bash virus_assembly_bwa.sh -p test -r your_path/viral.genomic.fna -o test/ -t 8
   for generating a report:
   
        bash report.sh -p test -o test/ -m default
        
  


# Usage
pipeline.sh

requests:
  1. -1: input .fastq file or .fastq.gz file for paired-end reads sequence 1
  2. -2: input .fastq file or .fastq.gz file for paired-end reads sequence 2
  3. -h: host genome hisat2 index
  4. -l: the virus list 
  5. -k: kallisto index
  6. -r: reference viral genome

optional:
  1. -t: number of threads used; default 1
  2. -o: directory for output files; default is current path
  3. -m: assembly method used: assembly or assembly_kallisto or assembly_bwa; default assembly
  4. -i: path to store single viral genome reference; default will auto create an file refs in working path

preprocess.sh

request:
  1. -1: input .fastq file or .fastq.gz file for paired-end reads sequence 1
  2. -2: input .fastq file or .fastq.gz file for paired-end reads sequence 2
  3. -h: host genome hisat2 index

optional:
  1. -t: number of threads used; default 1
  2. -o: directory for output files; default is current path

virus_identification.sh

request:
  1. -p: prefix of input .fq files in preprocess output
  2. -l: the virus list 
  3. -k: kallisto index

optional:
  1. -t: number of threads used; default 1
  2. -o: directory for output files; default is current path

virus_assembly.sh; virus_assembly_kallisto.sh; virus_assembly_bwa.sh:

request:
  1. -p: prefix of input .fq files in preprocess output
  2. -r: reference viral genome
  3. -i: path to store single viral genome reference; default will auto create an file refs in working path

optional:
  1. -o: directory for output files; default is current path
  2. -t: number of threads used; default 1
  
report.sh:

request:
  1. -p: prefix of input .fq files in preprocess output
  2. -m: assembly method used: default or kallisto or bwa

optional:
  1. -o: directory for output files; default is current path




