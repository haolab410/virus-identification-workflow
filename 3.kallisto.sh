thread_num=1
#mkdir raw
#cat SRR_Acc_List.txt | xargs -I {} -P ${thread_num} sh -c "prefetch {}"
#mkdir qc
#mkdir trim
#mkdir pangolin/align
#cat SRR_Acc_List.txt | xargs -I {} -P ${thread_num} sh -c "mv /home/shiyangs/ncbi/public/sra/{}.sra /home/shiyangs/virus/raw/{}.sra"
#cat SRR_Acc_List.txt | xargs -I {} -P ${thread_num} sh -c "fastq-dump --gzip --split-3 -O ~/virus/raw/ -A /home/shiyangs/virus/raw/{}.sra; mv ~/virus/raw/_home_shiyangs_virus_raw_{}.sra_1.fastq.gz ~/virus/raw/{}.1.fq.gz; mv ~/virus/raw/_home_shiyangs_virus_raw_{}.sra_2.fastq.gz ~/virus/raw/{}.2.fq.gz"
cat SRR_Acc_List.human.txt | xargs -I {} -P ${thread_num} sh -c "mkdir ~/virus/human/{}.kallisto; bash ~/virus/pipeline.kallisto.sh ~/virus/raw/{}.1.fq.gz ~/virus/raw/{}.2.fq.gz ~/virus/human/{}.kallisto/"

#cat SRR_Acc_List.txt | xargs -I {} -P ${thread_num} sh -c "mkdir /home/shiyangs/virus/human/{}; fastqc -o /home/shiyangs/virus/human -t 20 ~/virus/raw/{}.1.fq.gz ~/virus/raw/{}.2.fq.gz; /public211/tingm/software/FxTools/bin/FxTools_Linux Fqtools fqcheck -i ~/virus/raw/{}.1.fq.gz ~/virus/raw/{}.2.fq.gz -o ~/virus/human/{}/out1 ~/virus/human/{}/out2; trim_galore --fastqc --paired ~/virus/raw/{}.1.fq.gz ~/virus/raw/{}.2.fq.gz -o ~/virus/; mkdir ~/virus/align/{}; "
#mkdir pangolin/FVE
#cat pangolin.txt | xargs -I {} -P ${thread_num} sh -c "mkdir ~/virus/pangolin/align/{}; hisat2 -x /home/shiyangs/virus/virus_hisat_ref/ref -q -1 ~/virus/pangolin/{}.1_val_1.fq.gz -2 ~/virus/pangolin/{}.2_val_2.fq.gz -p 8 -S ~/virus/pangolin/align/{}/human_mapped.sam 2> ~/virus/pangolin/align/{}/human.summary.txt; samtools view -f 4 -bS ~/virus/pangolin/align/{}/human_mapped.sam > ~/virus/pangolin/align/{}/human_mapped.bam; samtools sort ~/virus/pangolin/align/{}/human_mapped.bam -o ~/virus/pangolin/align/{}/human_unmapped.sort.bam"
#cat pangolin.txt | xargs -I {} -P ${thread_num} sh -c "java -jar ~/software/picard/picard.jar MarkDuplicates M=/home/shiyangs/virus/pangolin/align/{}/dupstats REMOVE_DUPLICATES=TRUE I=/home/shiyangs/virus/pangolin/align/{}/human_unmapped.sort.bam O=/home/shiyangs/virus/align/{}/human_unmapped_nodup.bam; /home/software/bam2fastq-1.1.0/bam2fastq ~/virus/pangolin/align/{}/human_unmapped_nodup.bam -o ~/virus/pangolin/align/{}/human_unmapped#.fq"
#mkdir pangolin/bwa_virus
#cat pangolin.txt | xargs -I {} -P ${thread_num} sh -c "java -cp /home/shiyangs/software/FastViromeExplorer-master/bin/ FastViromeExplorer -1 ~/virus/pangolin/align/{}/human_unmapped_1.fq -2 ~/virus/pangolin/align/{}/human_unmapped_2.fq -l /home/yhliu/viral_genomes/virus_obtaining/ncbi-viruses-list.txt -i /home/yhliu/viral_genomes/virus_obtaining/ncbi-virus-kallisto-index-k31.idx -o /home/shiyangs/virus/pangolin/FVE/{}"
#bash bwa_virus.sh
