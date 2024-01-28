#!/bin/bash
                             
usage() {                                      # Function: Print a help message.
  echo "Usage: $0 [ -h help] [ -i input ] [ -t treat ] [ -c quickstart_config.yml ]" 1>&2
}
exit_abnormal() {                              # Function: Exit with error.
  usage
  exit 1
}


while getopts ":h:i:t:c:" options; do              
                                                                                           
  case "${options}" in
    h|:)
      usage
      Help
      exit 0
      ;;
    i)                                         
      input=${OPTARG}
      if ! [[ -n $input ]] ; then                      
        echo "You didn't set the input sample"                                
      fi                          
      ;;
    t)                                         
      treat=${OPTARG}
      if ! [[ -n $treat ]] ; then                      
        echo "You didn't set the treat sample"                                
      fi                          
      ;;
    c)                                         
      config=${OPTARG}
      if ! [[ -n $config ]] ; then                      
        echo "You didn't set the config yml file"                                
      fi                          
      ;;
    \?) # incorrect option
      echo "Error: -${OPTARG} Invalid option"
      exit_abnormal
      ;;
  esac
done

shift $(($OPTIND - 1))

quickstart_config=$(cat $config) 
quickstart_array=($(echo $quickstart_config | tr ":" "\n"))

echo -e "==Starting a quickstart job...==:\n"

# group: "Mammal"#1
# genome: "Homo_sapiens"#3
# assembly: "hg38"#5
# se_input: "/public/home/chenzr/PSI_Seq_brainCell/A1-A12-totalRNA-result/A1.cutadapt.extendedFrags.collapse.cutBarcodes.fa.gz"#7
# se_treat: "/public/home/chenzr/PSI_Seq_brainCell/A1-A12-totalRNA-result/A2.cutadapt.extendedFrags.collapse.cutBarcodes.fa.gz"#9
# se_input_prefix: "A1"#11
# se_treat_prefix: "A2"#13
# se_br5_input: "0"#15
# se_br3_input: "0"#17
# se_ad5_input: ""#19
# se_ad3_input: ""#21
# se_br5_treat: "0"#23
# se_br3_treat: "0"#25
# se_ad5_treat: ""#27
# se_ad3_treat: ""#29
# sites_identification: "true"#31
# sites_annotation: "true"#33
# sites_target_prediction: "true"#35
# br_remove_input: "false"#37
# br_remove_treat: "false"#39
# ad_remove_input: "false"#41
# ad_remove_treat: "false"#43
# svm: "true"#45
# metagene: "false"#47
# ann: "false"#49
# user_defined: "false"#51
# treatpre_Fold_thres: "0"#53
# preFold_FC_thres: "0"#55
# treataft_Fold_thres: "0"#57
# aftFold_FC_thres: "0"#59
# treatstoprate_thres: "0"#61
# stoprate_FC_thres: "0"#63

#cutadapt-se
echo -e "ad_remove_input:${quickstart_array[41]}"
if [[ "${quickstart_array[41]}" =~ "true"  ]] ;
then
    echo -e "Getting input trimmed fastq...\n"
    snakemake -s snakemake/Snakefile --cores 8 snakemake/trimmed/input/"${input}".fastq
else
    echo -e "Moving input trimmed fastq...\n"
    mv snakemake/reads/input/"${input}".fastq snakemake/trimmed/input
fi

echo -e "ad_remove_treat:${quickstart_array[43]}"
if [[ "${quickstart_array[43]}" =~ "true"  ]] ;
then
    echo -e "Getting treat trimmed fastq...\n"
    snakemake -s snakemake/Snakefile --cores 8 snakemake/trimmed/treat/"${treat}".fastq
else
    echo -e "Moving treat trimmed fastq...\n"
    mv snakemake/reads/treat/"${treat}".fastq snakemake/trimmed/treat
fi

#STAR-index
temp=$(echo -e "${quickstart_array[5]}" | sed "s/\"//g"  )
DIRECTORY="snakemake/genome/${temp}"
if [ -d "$DIRECTORY" ]; then
  echo "$DIRECTORY is not empty, STAR index have been built"
  
else
  echo -e "No STAR index in $DIRECTORY, generating STAR index...\n"
  snakemake -s snakemake/Snakefile --cores 8 snakemake/genome/"${temp}"
fi

#rtsSeeker result
echo -e "sites_identification:${quickstart_array[31]}"
if [[ "${quickstart_array[31]}" =~ "true"  ]] ;
then
    echo -e "Getting rtsSeeker result for Ψ sites identification...\n"
    snakemake -s snakemake/Snakefile --cores 8 snakemake/output/sites_identification/"${input}"_versus_"${treat}"/"${input}"_versus_"${treat}".bed
fi

#rtsSeeker ann
echo -e "ann:${quickstart_array[49]}"
if [[ "${quickstart_array[49]}" =~ "true"  ]] ;
then
    echo -e "Getting ann result for Ψ sites identification...\n"
    snakemake -s snakemake/Snakefile --cores 8 snakemake/output/ann/"${input}"_versus_"${treat}"/"${input}"_versus_"${treat}"_ann_psi_prediction.bed
fi

#rtsSeeker svm
echo -e "svm:${quickstart_array[45]}"
if [[ "${quickstart_array[45]}" =~ "true"  ]] ;
then
    echo -e "Getting svm result for Ψ sites identification...\n"
    snakemake -s snakemake/Snakefile --cores 8 snakemake/output/svm/"${input}"_versus_"${treat}"/"${input}"_versus_"${treat}"_svm_psi_prediction.bed
fi

#rtsSeeker user-defined
echo -e "user-defined:${quickstart_array[51]}"
if [[ "${quickstart_array[51]}" =~ "true"  ]] ;
then
    echo -e "Getting user-defined result for Ψ sites identification...\n"
    snakemake -s snakemake/Snakefile --cores 8 snakemake/output/user_defined/"${input}"_versus_"${treat}"/"${input}"_versus_"${treat}"_user_defined_psi_prediction.bed
fi

#bedAnnotator
echo -e "sites_annotation:${quickstart_array[33]}"
if [[ "${quickstart_array[33]}" =~ "true"  ]] ;
then
    echo -e "Getting bedAnnotator result for Ψ sites annotation...\n"
    snakemake -s snakemake/Snakefile --cores 8 snakemake/output/sites_annotation/"${input}"_versus_"${treat}"/"${input}"_versus_"${treat}"_add_seq_group_uniq.bed
fi

#metagene
echo -e "metagene:${quickstart_array[47]}"
if [[ "${quickstart_array[47]}" =~ "true"  ]] ;
then
    echo -e "Getting metagene result for Ψ sites annotation...\n"
    snakemake -s snakemake/Snakefile --cores 8 snakemake/output/sites_annotation/"${input}"_versus_"${treat}"/"${input}"_versus_"${treat}"_metagene_pseudoU_norm_length.pdf
fi

#ACAscan
echo -e "sites_target_prediction:${quickstart_array[35]}"
if [[ "${quickstart_array[35]}" =~ "true"  ]] ;
then
    echo -e "Getting ACAscan result for Ψ sites snoRNA target prediction...\n"
    snakemake -s snakemake/Snakefile --cores 8 snakemake/output/sites_target_prediction/"${input}"_versus_"${treat}"/"${input}"_versus_"${treat}"_fa.out
    echo -e "Getting PUSscan result for Ψ sites PUS target prediction...\n"
    snakemake -s snakemake/Snakefile --cores 8 snakemake/output/sites_target_prediction/"${input}"_versus_"${treat}"/overall_multinomialnb_model_prediction.txt
fi

#emit done signal
echo -e "==All done!=="
