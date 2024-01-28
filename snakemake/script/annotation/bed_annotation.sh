#! /bin/bash

usage() {
    echo 'Usage: ./'$0 ' -i bedfile -a bed6file -b ann_bool -c svm_bool -d user_defined_bool -e output_path'  
    exit -1
}

while getopts 'i:a:b:c:d:e:' OPT; do
    case $OPT in
        i) bedfile="$OPTARG";;
        a) bed6file="$OPTARG";;
        b) ann="$OPTARG";;
        c) svm="$OPTARG";;
        d) user_defined="$OPTARG";;
        e) output_path="$OPTARG";;
        h) usage;;
        ?) usage;;
    esac
done

output_path=${output_path%_add_seq_group_uniq.bed}

if [[ "$ann" == "true" ]]
then #ANN
    echo -e "ann mode is selected!"

    echo "generating ${output_path}_anno.bed"
    cat $bedfile > ${output_path}.psi.bed
    cut -f 1-6 $bedfile > ${output_path}.bed6
    bedAnnotator_cmd1="./snakemake/script/annotation/bedAnnotator -s 1 --anno $bed6file --bed ${output_path}.bed6 -o ${output_path}_anno.bed"
    bedAnnotator_cmd2="./snakemake/script/annotation/bedAnnotator -s 1 --anno $bed6file --bed $bedfile -o ${output_path}_anno_append.bed"
    eval $bedAnnotator_cmd1 &
    eval $bedAnnotator_cmd2 &

    wait

    #add gene biotype
    perl $(dirname "$0")/add_biotype.pl ${output_path}_anno.bed

    #add sequence and remove redundancy
    Rscript $(dirname "$0")/bed_annotation.r -f ${output_path}_anno.biotype.bed -g $bedfile -s $(dirname "$0")/hg38.psiU.SingleSites.bed -e $(dirname "$0")/human.hg38.Pseudo.result.col29.xlsx -i ./snakemake/genome/hg38.fa -j $(dirname "$0")/snoRNAbase_rRNA.fa -k $(dirname "$0")/snoRNAbase_snRNA.fa -o ${output_path}

    echo "generate ${output_path}_anno_stopRatioFC.bed"
    sort -r -g -k5 -k7 ${output_path}_anno.bed > ${output_path}_anno_stopRatioFC.bed
    awk 'BEGIN{print "#chrom\tchromStart\tchromEnd\tname\tfoldChange\tstrand\tgeneName\tgeneType\tolBaseNum\tqueryCov\tsampCov\tupDist\tdownDist"}1' ${output_path}_anno_stopRatioFC.bed > ${output_path}_anno_stopRatioFC.txt

    echo "generate ${output_path}_uniq_chrloc.txt (Venn diagram input)"
    cat ${output_path}_add_seq_group_uniq.bed |awk 'FS=OFS="\t" {print ""$1"_"$2"_"$3"_"$6"_"$7""}' > ${output_path}_uniq_chrloc.txt

    echo "generate ${output_path}_uniq_chrloc.txt (Venn diagram input)"
    cat ${output_path}_add_seq_group_uniq.bed |awk 'FS=OFS="\t" {print $52}' > ${output_path}_uniq_seq.txt

    echo -e "Finished: bedAnnotator done\n" 

    echo -e "bedAnnotator result in $(dirname ${output_path})"

    mupdf-x11 ${output_path}_gene_biotype_piechart.pdf &> /dev/null

elif [[ "$svm" == "true" ]]
then #SVM
    echo -e "SVM mode is selected!"

    echo "generating ${output_path}_anno.bed" 
    cat $bedfile > ${output_path}.psi.bed
    cut -f 1-6 $bedfile > ${output_path}.bed6
    bedAnnotator_cmd1="./snakemake/script/annotation/bedAnnotator -s 1 --anno $bed6file --bed ${output_path}.bed6 -o ${output_path}_anno.bed"
    bedAnnotator_cmd2="./snakemake/script/annotation/bedAnnotator -s 1 --anno $bed6file --bed $bedfile -o ${output_path}_anno_append.bed"
    eval $bedAnnotator_cmd1 &
    eval $bedAnnotator_cmd2 &

    wait

    #add gene biotype
    perl $(dirname "$0")/add_biotype.pl ${output_path}_anno.bed

    #add sequence and remove redundancy
    Rscript $(dirname "$0")/bed_annotation.r -f ${output_path}_anno.biotype.bed -g $bedfile -s $(dirname "$0")/hg38.psiU.SingleSites.bed -e $(dirname "$0")/human.hg38.Pseudo.result.col29.xlsx -i ./snakemake/genome/hg38.fa -j $(dirname "$0")/snoRNAbase_rRNA.fa -k $(dirname "$0")/snoRNAbase_snRNA.fa -o ${output_path}

    echo "generate ${output_path}_anno_stopRatioFC.bed"
    sort -r -g -k5 -k7 ${output_path}_anno.bed > ${output_path}_anno_stopRatioFC.bed
    awk 'BEGIN{print "#chrom\tchromStart\tchromEnd\tname\tfoldChange\tstrand\tgeneName\tgeneType\tolBaseNum\tqueryCov\tsampCov\tupDist\tdownDist"}1' ${output_path}_anno_stopRatioFC.bed > ${output_path}_anno_stopRatioFC.txt

    echo "generate ${output_path}_uniq_chrloc.txt (Venn diagram input)"
    cat ${output_path}_add_seq_group_uniq.bed |awk 'FS=OFS="\t" {print ""$1"_"$2"_"$3"_"$6"_"$7""}' > ${output_path}_uniq_chrloc.txt

    echo "generate ${output_path}_uniq_chrloc.txt (Venn diagram input)"
    cat ${output_path}_add_seq_group_uniq.bed |awk 'FS=OFS="\t" {print $52}' > ${output_path}_uniq_seq.txt

    echo -e "Finished: bedAnnotator done\n" 

    echo -e "bedAnnotator result in $(dirname ${output_path})"

    mupdf-x11 ${output_path}_gene_biotype_piechart.pdf &> /dev/null

elif [[ "$user_defined" == "true" ]]
then #User-defined
    echo -e "User-defined mode is selected!"

    echo "generating ${output_path}_anno.bed" 
    cat $bedfile > ${output_path}.psi.bed
    cut -f 1-6 $bedfile > ${output_path}.bed6
    bedAnnotator_cmd1="./snakemake/script/annotation/bedAnnotator -s 1 --anno $bed6file --bed ${output_path}.bed6 -o ${output_path}_anno.bed"
    bedAnnotator_cmd2="./snakemake/script/annotation/bedAnnotator -s 1 --anno $bed6file --bed $bedfile -o ${output_path}_anno_append.bed"
    eval $bedAnnotator_cmd1 &
    eval $bedAnnotator_cmd2 &

    wait

    #add gene biotype
    perl $(dirname "$0")/add_biotype.pl ${output_path}_anno.bed

    #add sequence and remove redundancy
    Rscript $(dirname "$0")/bed_annotation.r -f ${output_path}_anno.biotype.bed -g $bedfile -s $(dirname "$0")/hg38.psiU.SingleSites.bed -e $(dirname "$0")/human.hg38.Pseudo.result.col29.xlsx -i ./snakemake/genome/hg38.fa -j $(dirname "$0")/snoRNAbase_rRNA.fa -k $(dirname "$0")/snoRNAbase_snRNA.fa -o ${output_path}

    echo "generate ${output_path}_anno_stopRatioFC.bed"
    sort -r -g -k5 -k7 ${output_path}_anno.bed > ${output_path}_anno_stopRatioFC.bed
    awk 'BEGIN{print "#chrom\tchromStart\tchromEnd\tname\tfoldChange\tstrand\tgeneName\tgeneType\tolBaseNum\tqueryCov\tsampCov\tupDist\tdownDist"}1' ${output_path}_anno_stopRatioFC.bed > ${output_path}_anno_stopRatioFC.txt

    echo "generate ${output_path}_uniq_chrloc.txt (Venn diagram input)"
    cat ${output_path}_add_seq_group_uniq.bed |awk 'FS=OFS="\t" {print ""$1"_"$2"_"$3"_"$6"_"$7""}' > ${output_path}_uniq_chrloc.txt

    echo "generate ${output_path}_uniq_chrloc.txt (Venn diagram input)"
    cat ${output_path}_add_seq_group_uniq.bed |awk 'FS=OFS="\t" {print $52}' > ${output_path}_uniq_seq.txt

    echo -e "Finished: bedAnnotator done\n" 

    echo -e "bedAnnotator result in $(dirname ${output_path})"

    mupdf-x11 ${output_path}_gene_biotype_piechart.pdf &> /dev/null
fi

