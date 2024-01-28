import os
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

fq1 = snakemake.input.get("fq1")
assert fq1 is not None, "input-> fq1 is a required input parameter"
fq1 = (
    [snakemake.input.fq1]
    if isinstance(snakemake.input.fq1, str)
    else snakemake.input.fq1
)
#fq2 = snakemake.input.get("fq2")
#if fq2:
#    fq2 = (
#        [snakemake.input.fq2]
#        if isinstance(snakemake.input.fq2, str)
#        else snakemake.input.fq2
#    )
#    assert len(fq1) == len(
#        fq2
#    ), "input-> equal number of files required for fq1 and fq2"
input_str_fq1 = ",".join(fq1)
#input_str_fq2 = ",".join(fq2) if fq2 is not None else ""
#input_str = " ".join([input_str_fq1, input_str_fq2])

if fq1[0].endswith(".gz"):
    readcmd = "--readFilesCommand zcat"
else:
    readcmd = ""

outprefix = os.path.dirname(snakemake.output[0]) + "/"
#outprefix = snakemake.output[0]

shell(
    "STAR "
    "{extra} "
    "--runThreadN {snakemake.threads} "
    "--genomeDir {snakemake.params.index} "
    "--readFilesIn {input_str_fq1} "
    "{readcmd} "
    "--outFileNamePrefix {outprefix} "
    "--outStd Log "
    "--outSAMtype BAM SortedByCoordinate "
    "--alignEndsType EndToEnd "
    "--genomeLoad NoSharedMemory "
    "--limitBAMsortRAM 60000000000 "
    "--outFilterType BySJout "
    "--outFilterMultimapScoreRange 0 "
    "--outFilterMultimapNmax 30 "
    "--outFilterMismatchNmax 15 "
    "--outFilterMismatchNoverLmax 0.1 "
    "--outFilterScoreMin 0 "
    "--outFilterScoreMinOverLread 0 "
    "--outFilterMatchNmin 15 "
    "--outFilterMatchNminOverLread 0.8 "
    "--alignIntronMin 20 "
    "--alignIntronMax 1000000 "
    "--alignMatesGapMax 1000000 "
    "--seedSearchStartLmax 15 "
    "--seedSearchStartLmaxOverLread 1 "
    "--seedSearchLmax 0 "
    "--seedMultimapNmax 20000 "
    "--seedPerReadNmax 1000 "
    "--seedPerWindowNmax 100 "
    "--seedNoneLociPerWindow 20 "
    "--alignSJoverhangMin 8 "
    "--alignSJDBoverhangMin 1 "
    "--quantMode TranscriptomeSAM GeneCounts "
    "--outSAMmode Full "
    "--outSAMattributes All "
    "--outSAMunmapped None "
    "--outSAMorder Paired "
    "--outSAMprimaryFlag AllBestScore "
    "--outSAMreadID Standard "
    "--outReadsUnmapped Fastx "
    "--alignEndsProtrude 150 ConcordantPair "
    "--limitOutSJcollapsed 5000000 "
    "{log}"
)
