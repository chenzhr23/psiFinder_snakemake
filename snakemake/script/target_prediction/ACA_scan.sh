#! /bin/bash

if [[ $# -ne 3 ]]
then
    echo 'Usage: ./'$0 ' bedfile acaboxseq output_path'  
    exit 1
fi 

bedfile=$1
acaboxseq=$2
output_path=$3
output_path=${output_path%_fa.out}

cat ${bedfile%_add_seq_group_uniq.bed}_anno.bed | cut -f 1-6 > ${bedfile%_add_seq_group_uniq.bed}_anno.bed6
cat ${bedfile%_add_seq_group_uniq.bed}_anno.bed |awk 'FS=OFS="\t" {print ""$1"_"$2"_"$3"_"$6"",$5,$7}' > ${bedfile%_add_seq_group_uniq.bed}_anno.bed6.append
bedtools intersect -a ${bedfile%_add_seq_group_uniq.bed}.psi.bed -b ${bedfile%_add_seq_group_uniq.bed}_anno.bed6 -s > ${output_path}_anno_info.bed
bedtools intersect -a ${output_path}_anno_info.bed -b $(dirname "$0")/hg38.psiU.SingleSites.bed -wb|awk 'FS=OFS="\t" {print ""$1"_"$2"_"$3"_"$6"",$43,$39}' > ${output_path}.known.target
bedtools intersect -a ${output_path}_anno_info.bed -b $(dirname "$0")/hg38.psiU.SingleSites.bed -wb -v |awk 'FS=OFS="\t" {print ""$1"_"$2"_"$3"_"$6"","unknown",$39}' > ${output_path}.unknown.target
cat ${output_path}.known.target ${output_path}.unknown.target > ${output_path}_pseudoU.sites
sed -i -e 's/T/U/g' ${output_path}_pseudoU.sites

wait

echo "generating ACAscan result(_fa.out)" 
$(dirname "$0")/ACAscan -f $acaboxseq -m ${output_path}_pseudoU.sites >${output_path}_fa.out

echo "generating text information for ACAscan result"
Rscript $(dirname "$0")/ACA_scan.r -f ${output_path}_fa.out -a ${bedfile%_add_seq_group_uniq.bed}_anno.bed6.append -b $(dirname "$0")/human_hg38_snoRNABase_snoDB_rmRepeat_addorphaninfo.csv -c $(dirname "$0")/hg38.psiU.SingleSites.bed -o ${output_path}

# mupdf-x11 ${output_path}_pseudoU_Orphan_piechart.pdf &> /dev/null
# mupdf-x11 ${output_path}_pseudoU_snoRNA_guided_RNA_piechart.pdf &> /dev/null

echo -e "Finished: ACAscan done!\n" 
echo -e "ACAscan result in $(dirname ${output_path})"
