from snakemake.shell import shell
import re


n = len(snakemake.input)
assert n == 1, "Input must contain 1 (single-end) element."

extra = snakemake.params.get("extra", "")
adapters = snakemake.params.get("adapters", "")
se_ad5_input = snakemake.params.get("se_ad5_input", "")
se_ad3_input = snakemake.params.get("se_ad3_input", "")
se_ad5_treat = snakemake.params.get("se_ad5_treat", "")
se_ad3_treat = snakemake.params.get("se_ad3_treat", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

assert (
    extra != "" or adapters != ""
), "No options provided to cutadapt. Please use 'params: adapters=' or 'params: extra='."

pattern = re.compile('[agctnAGCTN]+')
# match1 = pattern.findall(se_ad5_input)
# match2 = pattern.findall(se_ad3_input)
# match3 = pattern.findall(se_ad5_treat)
# match4 = pattern.findall(se_ad3_treat)
match1 = bool(re.match(pattern,se_ad5_input))
match2 = bool(re.match(pattern,se_ad3_input))
match3 = bool(re.match(pattern,se_ad5_treat))
match4 = bool(re.match(pattern,se_ad3_treat))

if match1 == True and match2 == True:        
    print('start running cutadapt process for input sample...')
    shell(
    "cutadapt"
    " {snakemake.params.adapters}"
    " {snakemake.params.extra}"
    " -j {snakemake.threads}"
    " -o {snakemake.output.fastq}"
    " {snakemake.input[0]}"
    " > {snakemake.output.qc} {log}")    
elif match3 == True and match4 == True:
    print('start running cutadapt process for treat sample...')
    shell(
    "cutadapt"
    " {snakemake.params.adapters}"
    " {snakemake.params.extra}"
    " -j {snakemake.threads}"
    " -o {snakemake.output.fastq}"
    " {snakemake.input[0]}"
    " > {snakemake.output.qc} {log}")
elif match1 == True and match2 == False:
    print('start running cutadapt process for input sample on 5_end adapter...')
    shell(
    "cutadapt"
    " -a {snakemake.params.se_ad5_input}"
    " {snakemake.params.extra}"
    " -j {snakemake.threads}"
    " -o {snakemake.output.fastq}"
    " {snakemake.input[0]}"
    " > {snakemake.output.qc} {log}")
elif match1 == False and match2 == True:
    print('start running cutadapt process for input sample on 3_end adapter...')
    shell(
    "cutadapt"
    " -g {snakemake.params.se_ad3_input}"
    " {snakemake.params.extra}"
    " -j {snakemake.threads}"
    " -o {snakemake.output.fastq}"
    " {snakemake.input[0]}"
    " > {snakemake.output.qc} {log}")
elif match3 == True and match4 == False:
    print('start running cutadapt process for treat sample on 5_end adapter...')
    shell(
    "cutadapt"
    " -a {snakemake.params.se_ad5_treat}"
    " {snakemake.params.extra}"
    " -j {snakemake.threads}"
    " -o {snakemake.output.fastq}"
    " {snakemake.input[0]}"
    " > {snakemake.output.qc} {log}")
elif match3 == False and match4 == True:
    print('start running cutadapt process for treat sample on 3_end adapter...')
    shell(
    "cutadapt"
    " -g {snakemake.params.se_ad3_treat}"
    " {snakemake.params.extra}"
    " -j {snakemake.threads}"
    " -o {snakemake.output.fastq}"
    " {snakemake.input[0]}"
    " > {snakemake.output.qc} {log}")
# else:
#     print('no adapters are provided, run default cutadapt process with barcode...')
#     shell(
#     "cutadapt"
#     " {snakemake.params.extra}"
#     " -j {snakemake.threads}"
#     " -o {snakemake.output.fastq}"
#     " {snakemake.input[0]}"
#     " > {snakemake.output.qc} {log}")            

# else:
#     print('no adapters are provided, skipping cutadapt process...')
#     shell(
#     "mv"
#     " {snakemake.input[0]}"	
#     " {snakemake.output.fastq}")

#     shell(
#     "echo"
#     " skip {snakemake.input[0]} adapter trimming because no adapters are provided!"	
#     " > {snakemake.output.qc}")                 

