outdir=.
while getopts ":i:o:" opt
do
	  case $opt in
		    i)
        i="$OPTARG"
		    ;;
		    o)
        outdir="$OPTARG"
		    ;;
        ?)
        echo "help"
        echo "-i       lineages-XXXX-XX-XX.csv generate from NCBItax2lin"
        echo "-o       output directory"
        exit 1;;
    esac
done

echo ${i}
cut the cols needed for our purpose
awk -F ',' '{print $1"\t"$2";"$3";"$4";"$5";"$6";"$7";"$8}' ${i} |sort -k 1b,1 |awk -F '\t' '{print $1":"$2}' >${outdir}/lineages-short.csv

#download virus genome from ncbi
cd ${outdir}
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.1.genomic.fna.gz

gzip -d viral.*.1.genomic.fna.gz
cat viral.*.1.genomic.fna >viral.genomic.fna
rm viral.*.1.genomic.fna
samtools faidx ${outdir}/viral.genomic.fna


# downloading accession2taxid file
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
gzip -d nucl_gb.accession2taxid.gz
grep -E 'AC|NC' nucl_gb.accession2taxid >nucl_gb.accession2taxid_virus
awk '{print $1"\t"$3}' nucl_gb.accession2taxid_virus |sort -k 1b,1 >nucl_gb.accession2taxid_virus_short

# collecting gb.acc from virus bank
grep '>' ${outdir}/viral.genomic.2020.01.27.fna|awk '{print $1}'|sed 's/>//'|awk -F '.' '{print $1}' >temp.1
grep '>' ${outdir}/viral.genomic.2020.01.27.fna|awk '{print $1}'|sed 's/>//' >temp.2
paste -d '\t' temp.1 temp.2 |sort -k 1b,1 >GBaccs.list.txt

# join files
join -a 1 -o 2.2 -o 1.1 -o 1.2 GBaccs.list.txt nucl_gb.accession2taxid_virus_short |sort -k1b,1|awk '{print $1":"$2":"$3}' >GBaccToTaxID.list
join -a 1 -t ':' -o 1.2 -o 1.3 -o 1.1 -o 2.2 GBaccToTaxID.list ${outdir}/lineages-short.csv|sed 's/:/\t/g' >${outdir}/taxons_virus.txt


grep '>' viral.genomic.2020.01.27.fna|sed 's/>\([A-Z|0-9|.|_]*\)\s\(.*\)/\1\t\2/'|sort -k1 >${outdir}/temp.ids_names.txt

awk '{print $1"\t"$2}' ${outdir}/viral.genomic.2020.01.27.fna.fai|sort -k1 > ${outdir}/temp.ids_length.txt

awk -F '\t' '{print $2"\t"$4}' ${outdir}/taxons_virus.txt |sort -k1 >${outdir}/temp.ids_taxons.txt

paste -d '\t' ${outdir}/temp.ids_names.txt ${outdir}/temp.ids_taxons.txt ${outdir}/temp.ids_length.txt|awk -F '\t' '{print $1"\t"$2"\t"$4"\t"$6}' >${outdir}/viruses-list.txt

rm ${outdir}/temp.ids_*
kallisto index -i ncbi-virus-kallisto-index-k31.idx viral.genomic.fna