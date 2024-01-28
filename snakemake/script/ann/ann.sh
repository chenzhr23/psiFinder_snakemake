#! /bin/bash

if [[ $# -ne 2 ]]
then
    echo 'Usage: '$0 ' rtsSeeker.bed output_path'  
    exit 1
fi 

out=$1
output_path=$2
output_path=${output_path%_ann_psi_prediction.bed}

echo "rtsSeeker ann mode: start filtering rtsSeeker result..."
awk 'FS=OFS="\t" {if($1~/^chr[0-9|a-z|A-Z]*$/ && $10=="T" && $14>10){print $0}}' $out > ${output_path}_ann_filt_totalRNA.bed
bedtools intersect -a ${output_path}_ann_filt_totalRNA.bed -b $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -s > ${output_path}_knowpse.bed
bedtools intersect -a ${output_path}_ann_filt_totalRNA.bed -b $(dirname "$0")/rrna_chr21.bed -s > ${output_path}_rrna.bed
bedtools intersect -a ${output_path}_rrna.bed -b $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -s -v > ${output_path}_notpse.bed
cat ${output_path}_knowpse.bed |awk 'FS=OFS="\t" {print $0,"1"}' > ${output_path}_knowpse.txt
cat ${output_path}_notpse.bed |awk 'FS=OFS="\t" {print $0,"0"}' > ${output_path}_notpse.txt
echo "getting roc_plot.txt"
cat ${output_path}_knowpse.txt ${output_path}_notpse.txt > ${output_path}_ann_plot.txt

nohup Rscript $(dirname "$0")/ann_build_totalRNA.r -f ${output_path}_ann_plot.txt -k ${output_path}_ann_filt_totalRNA.bed -r $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -s $(dirname "$0")/hg38.psiU.SingleSites.bed  -o ${output_path} > ${output_path}_ann_evaluation_totalRNA.log 2>&1 &
wait

cat ${output_path}_ann_evaluation_totalRNA.log
# mupdf-x11 ${output_path}_six_variables_rRNA_violinplot.pdf &> /dev/null
mupdf-x11 ${output_path}_ann_roc_test_data_plot.pdf &> /dev/null
mupdf-x11 ${output_path}_ann_evaluation.pdf &> /dev/null

echo "total RNA: ann build and predict programs are done!"
echo -e "total RNA: ann prediction result in $(dirname ${output_path}_ann_psi_prediction.bed)"
