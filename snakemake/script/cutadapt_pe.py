from snakemake.shell import shell
import re


n = len(snakemake.input)
assert n == 2, "Input must contain 2 (paired-end) elements."

extra = snakemake.params.get("extra", "")
adapters = snakemake.params.get("adapters", "")
pe_ad5_read1 = snakemake.params.get("pe_ad5_read1", "")
pe_ad3_read1 = snakemake.params.get("pe_ad3_read1", "")
pe_ad5_read2 = snakemake.params.get("pe_ad5_read2", "")
pe_ad3_read2 = snakemake.params.get("pe_ad3_read2", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

assert (
    extra != "" or adapters != ""
), "No options provided to cutadapt. Please use 'params: adapters=' or 'params: extra='."

pattern = re.compile('[agctAGCT]+')
match1 = pattern.findall(pe_ad5_read1)
match2 = pattern.findall(pe_ad3_read1)
match3 = pattern.findall(pe_ad5_read2)
match4 = pattern.findall(pe_ad3_read2)

if match1 and match2 and match3 and match4:        
    print('start running cutadapt process...')
    shell(
    "cutadapt"
    " {snakemake.params.adapters}"
    " {snakemake.params.extra}"
    " -o {snakemake.output.fastq1}"
    " -p {snakemake.output.fastq2}"
    " -j {snakemake.threads}"
    " {snakemake.input}"
    " > {snakemake.output.qc} {log}")   
else:
    print('no adapters are provided, skipping cutadapt process...')
    shell(
    "cp"
    " {snakemake.input[0]}"	
    " {snakemake.output.fastq1}")

    shell(
    "cp"
    " {snakemake.input[1]}" 
    " {snakemake.output.fastq2}")

    shell(
    "echo"
    " skip {snakemake.input[0]} {snakemake.input[1]} adapter trimming because no adapters are provided!"	
    " > {snakemake.output.qc}")