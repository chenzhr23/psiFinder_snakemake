#! /bin/bash

if [[ $# -ne 5 ]]
then
    echo 'Usage: ./'$0 ' input_sample input_path treat_sample treat_path log_file'
    exit 1
fi

input_sample=$1
input_path=$2
treat_sample=$3
treat_path=$4
log_file=$5


echo -e "clear psiFinder.sh history"
echo -e "Your quickstart job is starting, please wait..." > $log_file

#.fastq.gz
if [[ "$input_sample" == *.fastq.gz ]] ;
then
    echo -e "$(tput setaf 2)input: unzipping and redirecting fastq file...\n$(tput sgr 0)" >> $log_file
    gunzip -c $input_sample  > $input_path &
fi

if [[ "$treat_sample" == *.fastq.gz ]] ;
then
    echo -e "$(tput setaf 2)treat: unzipping and redirecting fastq file...\n$(tput sgr 0)" >> $log_file
    gunzip -c $treat_sample  > $treat_path &
fi

#.fasta.gz
if [[ "$input_sample" == *.fasta.gz ]] ;
then
    echo -e "$(tput setaf 2)input: unzipping and redirecting fasta file...\n$(tput sgr 0)" >> $log_file
    gunzip -c $input_sample  > $input_path
    seqtk seq -F '#' $input_path > ${input_path%.fasta}.fastq
    rm $input_path
fi

if [[ "$treat_sample" == *.fasta.gz ]] ;
then
    echo -e "$(tput setaf 2)treat: unzipping and redirecting fasta file...\n$(tput sgr 0)" >> $log_file
    gunzip -c $treat_sample  > $treat_path
    seqtk seq -F '#' $treat_path > ${treat_path%.fasta}.fastq
    rm $treat_path
fi

#.fq.gz
if [[ "$input_sample" == *.fq.gz ]] ;
then
    echo -e "$(tput setaf 2)input: unzipping and redirecting fastq file...\n$(tput sgr 0)" >> $log_file
    gunzip -c $input_sample  > $input_path &
fi

if [[ "$treat_sample" == *.fq.gz ]] ;
then
    echo -e "$(tput setaf 2)treat: unzipping and redirecting fastq file...\n$(tput sgr 0)" >> $log_file
    gunzip -c $treat_sample  > $treat_path &
fi

#.fa.gz
if [[ "$input_sample" == *.fa.gz ]] ;
then
    echo -e "$(tput setaf 2)input: unzipping and redirecting fasta file...\n$(tput sgr 0)" >> $log_file
    gunzip -c $input_sample  > $input_path
    seqtk seq -F '#' $input_path > ${input_path%.fasta}.fastq
    rm $input_path
fi

if [[ "$treat_sample" == *.fa.gz ]] ;
then
    echo -e "$(tput setaf 2)treat: unzipping and redirecting fasta file...\n$(tput sgr 0)" >> $log_file
    gunzip -c $treat_sample  > $treat_path
    seqtk seq -F '#' $treat_path > ${treat_path%.fasta}.fastq
    rm $treat_path
fi

#.fastq
if [[ "$input_sample" == *.fastq ]] ;
then
    echo -e "$(tput setaf 2)input: redirecting fastq file...\n$(tput sgr 0)" >> $log_file
    cp -R $input_sample $input_path &
fi

if [[ "$treat_sample" == *.fastq ]] ;
then
    echo -e "$(tput setaf 2)treat: redirecting fastq file...\n$(tput sgr 0)" >> $log_file
    cp -R $treat_sample $treat_path &
fi

#.fasta
if [[ "$input_sample" == *.fasta ]] ;
then
    echo -e "$(tput setaf 2)input: redirecting fasta file...\n$(tput sgr 0)" >> $log_file
    cp -R $input_sample $input_path
    seqtk seq -F '#' $input_path > ${input_path%.fasta}.fastq
    rm $input_path
fi

if [[ "$treat_sample" == *.fasta ]] ;
then
    echo -e "$(tput setaf 2)treat: redirecting fasta file...\n$(tput sgr 0)" >> $log_file
    cp -R $treat_sample $treat_path
    seqtk seq -F '#' $treat_path > ${treat_path%.fasta}.fastq
    rm $treat_path
fi

#.fq
if [[ "$input_sample" == *.fq ]] ;
then
    echo -e "$(tput setaf 2)input: redirecting fastq file...\n$(tput sgr 0)" >> $log_file
    cp -R $input_sample $input_path &
fi

if [[ "$treat_sample" == *.fq ]] ;
then
    echo -e "$(tput setaf 2)treat: redirecting fastq file...\n$(tput sgr 0)" >> $log_file
    cp -R $treat_sample $treat_path &
fi

#.fa
if [[ "$input_sample" == *.fa ]] ;
then
    echo -e "$(tput setaf 2)input: redirecting fasta file...\n$(tput sgr 0)" >> $log_file
    cp -R $input_sample $input_path
    seqtk seq -F '#' $input_path > ${input_path%.fasta}.fastq
    rm $input_path
fi

if [[ "$treat_sample" == *.fa ]] ;
then
    echo -e "$(tput setaf 2)treat: redirecting fasta file...\n$(tput sgr 0)" >> $log_file
    cp -R $treat_sample $treat_path
    seqtk seq -F '#' $treat_path > ${treat_path%.fasta}.fastq
    rm $treat_path
fi

wait


echo -e "Starting snakemake workflow..."
#snakemake
#fastq
if [[ "$(basename ${input_path})" == *.fastq ]] && [[ "$(basename ${treat_path})" == *.fastq ]];
then
    echo "check info: Input is fastq files"  >> $log_file
    nohup bash snakemake/psiFinder.sh -i $(basename ${input_path%.fastq}) -t $(basename ${treat_path%.fastq}) -c quickstart_config.yml >> $log_file 2>&1 &
    wait
else
    echo "check info: Input is not fastq files" >> $log_file
fi

#fasta
if [[ "$(basename ${input_path})" == *.fasta ]] && [[ "$(basename ${treat_path})" == *.fasta ]] ;
then
    echo "check info: Input is fasta files" >> $log_file
    nohup bash snakemake/psiFinder.sh -i $(basename ${input_path%.fasta}) -t $(basename ${treat_path%.fasta}) -c quickstart_config.yml >> $log_file 2>&1 &
    wait
else
    echo "check info: Input is not fasta files" >> $log_file
fi

echo -e "Finished snakemake workflow!" >> $log_file
