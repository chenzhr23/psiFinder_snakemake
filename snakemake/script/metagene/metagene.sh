#! /bin/bash

usage() {
	echo 'Usage: ./'$0 ' -i bedfile -g genome -a gtf [option] '
    exit -1
}

bedfile=''
genome=''
gtf=''
annofile=''

while getopts 'i:g:a:' OPT; do
	case $OPT in
		i) bedfile="$OPTARG";;
		g) genome="$OPTARG";; 
		a) gtf="$OPTARG";;
		h) usage;;
		?) usage;;
    esac
done
 
#get data for metagene plot
echo -e "get metagene data"
if [ ! -f ${gtf%.gtf}_annot_sorted.bed ]
	then
	if [ ! -f "${gtf%.gtf}.genePred"  ]
	then
	echo "generate genePredfile"
	gtfToGenePred -genePredExt -geneNameAsName2 $gtf ${gtf%.gtf}.genePred
	else
	echo "genePredfile exist"
	fi
	echo "generate metagene annotation bed"
	perl ./snakemake/script/metagene/make_annot_bed.pl --genomeDir $(dirname ${gtf%.gtf})/genome/ --genePred ${gtf%.gtf}.genePred > ${gtf%.gtf}_annot.bed
	bedtools sort -i ${gtf%.gtf}_annot.bed > ${gtf%.gtf}_annot_sorted.bed
	rm ${gtf%.gtf}_annot.bed
	echo  "generate region sizes data"
	perl ./snakemake/script/metagene/size_of_cds_utrs.pl --annot ${gtf%.gtf}_annot_sorted.bed >${gtf%.gtf}_region_sizes.txt
fi
echo -e "metagene annotaion"
cut -f 1-6 ${bedfile%_add_seq_group_uniq.bed}_anno.bed > ${bedfile%_add_seq_group_uniq.bed}_anno.bed6
perl ./snakemake/script/metagene/annotate_bed_file.pl --bed ${bedfile%_add_seq_group_uniq.bed}_anno.bed6 --bed2 ${gtf%.gtf}_annot_sorted.bed >${bedfile%_add_seq_group_uniq.bed}.sorted_annot.bed
perl ./snakemake/script/metagene/rel_and_abs_dist_calc.pl --bed ${bedfile%_add_seq_group_uniq.bed}.sorted_annot.bed --regions ${gtf%.gtf}_region_sizes.txt >${bedfile%_add_seq_group_uniq.bed}_metagene.txt

echo -e "generate metagene plot from ${bedfile%_add_seq_group_uniq.bed}_metagene.txt"
Rscript ./snakemake/script/metagene/metagene.r -f ${bedfile%_add_seq_group_uniq.bed}_metagene.txt -o ${bedfile%_add_seq_group_uniq.bed}
cp ${bedfile%_add_seq_group_uniq.bed}_metagene_pseudoU_norm_length.pdf ./snakemake/metagene/$(basename ${bedfile%_add_seq_group_uniq.bed})
mupdf-x11 ${bedfile%_add_seq_group_uniq.bed}_metagene_pseudoU_norm_length.pdf &> /dev/null

echo "metagene plot end"
echo -e "metagene result in ./snakemake/$(basename ${bedfile%_add_seq_group_uniq.bed})"
