import sys,subprocess
packagename=sys.argv[1] + '/'
sample = sys.argv[2]
count = int(sys.argv[3])
flag = int(sys.argv[4])
import pandas as pd

sample = sample.strip('\n')
cov = pd.read_csv(packagename + "cov.summary.txt",sep="\t")
cov.head()
cov["Reads_hit"] = cov.apply(lambda x: x["Plus_reads"]+x["Minus_reads"], axis=1)
data1 = cov[["#ID","Avg_fold","Covered_percent","Reads_hit"]]
data1.columns = ["GI","Average_depth_of_coverage","%Coverage","Reads_hit"]
abundance = pd.read_csv(packagename + "FastViromeExplorer-final-sorted-abundance.tsv",sep="\t")
abundance.head()
abundance["genus"] = abundance.iloc[:,2].str.split(';').str[5]
abundance["species"] = abundance.iloc[:,2].str.split(';').str[6]
data2 = abundance[["#VirusIdentifier","species","genus"]]
data2.columns = ["GI","Species","genus"]
data = pd.merge(data1,data2,how = 'left',on='GI')
new_data = data[["Species","genus","GI","%Coverage","Reads_hit","Average_depth_of_coverage"]]
new_data.to_csv(packagename + "virus.distribution.csv",header=True,index=False)
    
cov = pd.read_csv(packagename + 'cov.summary.txt', sep='\t')
cov.head()
new_data = cov[["#ID","Length", "Read_GC"]]
new_data.to_csv(packagename + "stats.cons.csv", header=True, index=False)
    
    
outFile = open(packagename + 'qc.report.txt', 'w')
outFile.write('Samples\tReads_N\tBase_N\tAverage_Length\tGC(%)\tQ20(%)\tQ30(%)\n')
if(flag == 1):
    for k in range(0,count):
        j = str(int(k) + 1)    # the j's file of all input file when multiple input file need to be merged
        inputFile = open(packagename + j + '.out1.fqcheck')
        i = inputFile.readlines()
        o1 = i[0].strip('\n')
        o1 = o1.split(',')
        o1 = o1[1].split(' ')[1] + '\t' + o1[2].split(' ')[1] + '\t' + o1[4].split(':')[-1] + '\t'
        o2 = i[-1].split('\t')
        o2 = o2[1] + '\t' + o2[2] + '\t' + o2[3]
        outFile.write(j + '.' + sample + '_R1\t' + o1 + o2)
        inputFile = open(packagename + j + '.out2.fqcheck')
        i = inputFile.readlines()
        o1 = i[0].strip('\n')
        o1 = o1.split(',')
        o1 = o1[1].split(' ')[1] + '\t' + o1[2].split(' ')[1] + '\t' + o1[4].split(':')[-1] + '\t'
        o2 = i[-1].split('\t')
        o2 = o2[1] + '\t' + o2[2] + '\t' + o2[3]
        outFile.write(j + '.' + sample + '_R2\t' + o1 + o2)
else:
    for k in range(0,count):
        j = str(int(k) + 1)
        inputFile = open(packagename + j + '.out.fqcheck')
        i = inputFile.readlines()
        o1 = i[0].strip('\n')
        o1 = o1.split(',')
        o1 = o1[1].split(' ')[1] + '\t' + o1[2].split(' ')[1] + '\t' + o1[4].split(':')[-1] + '\t'
        o2 = i[-1].split('\t')
        o2 = o2[1] + '\t' + o2[2] + '\t' + o2[3]
        outFile.write(j + '.' + sample + '_R1\t' + o1 + o2)