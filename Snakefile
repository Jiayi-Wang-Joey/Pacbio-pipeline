import pandas as pd
from snakemake.utils import min_version
import json
import itertools
import glob
configfile: "config.yaml"
R = config["R"]


### SET UP ###
configfile: "/home/jiayiwang/pacbio/config.yaml"
min_version(config["snakemake_min_version"])
container: f"docker://condaforge/mambaforge:{config['mambaforge_version']}"
res_dir = "/home/jiayiwang/loggedfs_david_penton/"

### Wildcards ###
# SAMPLE = ["PB205_r0243-bcM0002", "PB205_r0243-bcM0003", "PB205_r0243-bcM0004",
#        "PB205_r0246_bcM0002", "PB205_r0246_bcM0003", "PB205_r0246_bcM0004",
#        "PB202_1_r0237_bcM0001", "PB202_2_r0238_bcM0002", "PB261_1_r0389_bcM0001", 
#        "PB214_2_r261_bcM0002", "PB214_1_r260_bcM0001"]
SAMPLE = ["PB261_2_r0392_bcM0002", "PB214_2_r261_bcM0002"]

MERGE = {"tumor1": ["PB202_1_r0237_bcM0001"],
         "normal1":["PB202_2_r0238_bcM0002"],
         "normal2": ["PB205_r0243-bcM0003", "PB205_r0246_bcM0003"],
         "normal3": ["PB205_r0243-bcM0002", "PB205_r0246_bcM0002"],
         "tumor3": ["PB205_r0246_bcM0004", "PB205_r0243-bcM0004", "PB258_3_r0362_bcM0004"],
         "normal4": ["PB214_2_r261_bcM0002"],
         "tumor4": ["PB214_1_r260_bcM0001"],
         "tumor6": ["PB261_1_r0389_bcM0001"],
         "normal6": ["PB261_2_r0392_bcM0002"]
         }

TISSUE = ["normal1","tumor1", "normal3", "tumor3", 
          "normal2", "normal4", "tumor4", "tumor6", "normal6"]

# TISSUE = ["normal4", "normal6"]


### RESULT ### 
# preprocessing
primer_removal =  expand(res_dir + "result/preprocessing/primer_removal/{sample}.fl.5p--3p.bam", sample=SAMPLE)
tag = expand(res_dir + "result/preprocessing/tags/{sample}.flt.bam", sample=SAMPLE)
refine = expand(res_dir + "result/preprocessing/refine/{sample}.fltnc.bam", sample=SAMPLE)
bcc = expand(res_dir + "result/preprocessing/BC_correction/{sample}.corrected.bam", sample=SAMPLE)
sort = expand(res_dir + "result/preprocessing/sorting/{sample}.sorted.bam",sample=SAMPLE)
dedup = expand(res_dir + "result/preprocessing/deduplication/{sample}.unmapped.bam", sample=SAMPLE)

# align to transcriptome
transcriptome = "data/genome/transcriptome.fa"
keep_bc = expand(res_dir + "result/align/transcriptome/keep_BC/{sample}.tags.bam", sample=SAMPLE)
fq = expand(res_dir + "result/align/transcriptome/bam2fq/{sample}.fastq.gz", sample=SAMPLE)
align = expand(res_dir + "result/align/transcriptome/run_minimap2_novel_transcriptome/{sample}.aligned.bam",sample=SAMPLE)
bc_tag = expand(res_dir + "result/align/transcriptome/add_BC_tags/{sample}.aligned.tagged.bam", sample=SAMPLE)
merge1 = expand(res_dir + "result/align/transcriptome/merge/{tissue}.aligned.tagged.bam", tissue=TISSUE)

# quantification
bc_sort = expand(res_dir+"result/quantify/transcriptome/sort_bams_by_barcode/{tissue}.aligned.sorted_barcode.bam", tissue=TISSUE)
oarfish = expand(res_dir+"result/quantify/transcriptome/oarfish/{tissue}.count.mtx", tissue=TISSUE)
sce = expand(res_dir+"result/quantify/transcriptome/gene_sce/{tissue}.rds", tissue=TISSUE)
#split = directory(expand(res_dir+"result/quantify/transcriptome/split_bams_by_barcode/{tissue}", tissue=TISSUE))
#barcode = expand(res_dir+"result/quantify/transcriptome/get_barcodes/{tissue}.barcodes.txt", tissue=TISSUE)
##collate_bam = expand(res_dir+"result/quantify/transcriptome/collate_bams/{tissue}.aligned.sorted_barcode.collate.bam", tissue=TISSUE)
#oarfish = directory(expand(res_dir+"result/quantify/transcriptome/run_oarfish/{tissue}", tissue=TISSUE))
transcript = expand(res_dir+"result/quantify/transcriptome/extract_transcript_gene_mapping/transcript_gene_mapping.tsv", tissue=TISSUE)
transcript_sce = expand(res_dir+"result/quantify/transcriptome/quantify_build_transcriptome_sce/{tissue}/transcriptome_sce.rds", tissue=TISSUE)
gene_sce = expand(res_dir+"result/quantify/transcriptome/quantify_build_genome_sce/{tissue}/genome_sce.rds", tissue=TISSUE)

# align to genome
align_genome = expand("/home/jiayiwang/loggedfs_david_penton/result/align/genome/run_minimap2_genome/{sample}.aligned.sorted.bam", sample=SAMPLE)
add_BC_tags_genome = expand("/home/jiayiwang/loggedfs_david_penton/result/align/genome/add_BC_tags/{merge}.aligned.tagged.bam", merge=MERGE)


preprocess = {
    "pmr": primer_removal,
    "tag": tag,
    "ref": refine,
    "bcc": bcc,
    "sort": sort,
    "dedup": dedup
}

transcriptome_align = {
    "tx": transcriptome,
    "keep_bc": keep_bc,
    "fq": fq,
    "align": align,
    "merge_tx": merge1
}

transcriptome_quantify = {
    # "bc_sort": bc_sort,
    # "oarfish": oarfish,
    "sce": sce
}

genome_res = {
    #"align_genome": align_genome,
    "add_BC_tags": add_BC_tags_genome
}




rule all:
    input:
        # [x for x in preprocess.values()],
        # [x for x in transcriptome_align.values()],
        [x for x in transcriptome_quantify.values()],
        #[x for x in genome_res.values()]
        #[x for x in variant_res.values()]

    #    expand(
    #        "results/quantify/transcriptome/quantify_build_genome_sce/{sample}/genome_sce.rds",
    #        sample=config["sample_names"],
    #    ),
    #    expand(
    #        "results/quantify/transcriptome/quantify_build_genome_sce/{sample}/transcriptome_sce.rds",
    #        sample=config["sample_names"],
    #    )


####################### preprocessing ########################
rule primer_removal:
    priority: 99
    input:
        reads="/home/jiayiwang/loggedfs_david_penton/data/kinnex/segmented/{sample}-segmented.bam",
        primers=config["primers"]
    output:
        res_dir + "result/preprocessing/primer_removal/{sample}.fl.5p--3p.bam"
    log:
        stderr="logs/preprocessing/primer_removal/{sample}.stderr"
    threads: config["lima_threads"]
    params:
        threads=config["lima_threads"],
        filename=res_dir + "result/preprocessing/primer_removal/{sample}.fl.bam"
    conda:
        "envs/preprocessing.yaml"
    shell:
        """
        lima {input.reads} {input.primers} {params.filename} --no-reports --num-threads {params.threads} --isoseq
        """


rule tag:
    priority: 98
    input:
        #"results/preprocessing/primer_removal/{sample}.fl.5p--3p.bam"
        rules.primer_removal.output
    output:
        res_dir + "result/preprocessing/tags/{sample}.flt.bam"
    log:
        stderr="logs/preprocessing/tags/{sample}.stderr"
    threads: config["tag_threads"]
    params:
        threads=config["tag_threads"],
    conda:
        "envs/preprocessing.yaml"
    shell:
        """
        isoseq tag {input} {output} --design T-12U-16B --num-threads {params.threads}
        """


rule refine:
    priority: 97
    input:
        reads=rules.tag.output,
        primers=config["primers"],
    output:
        res_dir + "result/preprocessing/refine/{sample}.fltnc.bam",
    log:
        stderr="logs/preprocessing/refine/{sample}.stderr",
    threads: config["refine_threads"]
    params:
        threads=config["refine_threads"],
    conda:
        "envs/preprocessing.yaml"
    shell:
        """
        isoseq refine {input.reads} {input.primers} {output} --require-polya --num-threads {params.threads}
        """


rule BC_correction:
    priority: 96
    input:
        reads=rules.refine.output,
        # barcodes=config["barcodes"],  # No need for this line now
    output:
        res_dir + "result/preprocessing/BC_correction/{sample}.corrected.bam",
    log:
        stderr="logs/preprocessing/BC_correction/{sample}.stderr",
    threads: config["bc_corr_threads"]
    params:
        threads=config["bc_corr_threads"],
        barcodes=lambda wildcards: "data/barcode/3M-3pgex-may-2023.REVCOMP.txt" 
        if wildcards.sample in ["PB214_1_r260_bcM0001", "PB214_2_r261_bcM0002", "PB261_1_r0389_bcM0001", "PB261_2_r0392_bcM0002"] 
        else "data/barcode/3M-february-2018-REVERSE-COMPLEMENTED.txt"
    conda:
        "envs/preprocessing.yaml"
    shell:
        """
        isoseq correct --barcodes {params.barcodes} \
            --num-threads {params.threads} {input.reads} {output}
        """



rule sort_bam:
    priority: 95
    input:
        rules.BC_correction.output,
    output:
        res_dir + "result/preprocessing/sorting/{sample}.sorted.bam",
    log:
        stderr="logs/preprocessing/sorting/{sample}.stderr",
    threads: config["sorting_threads"]
    params:
        threads=config["sorting_threads"],
        memory=config["sorting_memory"],
    conda:
        "envs/preprocessing.yaml"
    shell:
        """
        samtools sort -t CB -@ {params.threads} -m{params.memory}g {input} -o {output}
        """


# rule deduplication:
#     priority: 94
#     input:
#         rules.sort_bam.output
#     output:
#         res_dir + "result/preprocessing/deduplication/{sample}.unmapped.bam"
#     log:
#         stderr="logs/preprocessing/deduplication/{sample}.stderr"
#     threads: config["dedup_threads"]
#     params:
#         threads=config["dedup_threads"],
#     conda:
#         "envs/preprocessing.yaml"
#     shell:
#         """
#         isoseq groupdedup --num-threads {params.threads} {input} {output}
#         """

####################### alignment to transcriptome #######################
rule align_extract_transcriptome:
    input:
        # TODO: Replace with GENCODE GTF here
        transcriptome="data/genome/gencode.v47.chr_patch_hapl_scaff.annotation.gtf",
        genome="data/genome/GRCh38.p14.genome.fa"
    output:
        "data/genome/extract_transcriptome/transcriptome.fa"
    conda:
        "envs/gffread.yaml"
    shell:
        """
        gffread -w {output} -g {input.genome} {input.transcriptome}
        """


rule add_BC_tags_to_name:
    input:
        res_dir + "result/preprocessing/deduplication/{sample}.unmapped.bam"
    output:
        res_dir + "result/align/transcriptome/keep_BC/{sample}.tags.bam"
    log:
        stderr="logs/align/transcriptome/keep_BC/{sample}.stderr"
    conda:
        "envs/jvarkit.yaml"
    shell:
        """
        samjdk --samoutputformat BAM \
            -e '
            String cb = (String) record.getAttribute(\\"CB\\"); 
            String name = record.getReadName(); 
            record.setReadName(cb + \\"_\\" + name); 
            return record;' {input} > {output} 2> {log.stderr}
        """


rule bam2fq:
    input:
        res_dir + "result/align/transcriptome/keep_BC/{sample}.tags.bam"
    output:
        res_dir + "result/align/transcriptome/bam2fq/{sample}.fastq.gz"
    log:
        stderr="logs/align/transcriptome/bam2fq/{sample}.stderr"
    threads: config["convert_threads"]
    params:
        threads=config["convert_threads"]
    #conda:
    #    "envs/samtools.yaml"
    shell:
        """
        samtools bam2fq -@ {params.threads} {input} | gzip > {output}
        """


rule align_run_minimap2_novel_transcriptome:
    input:
        reads=rules.bam2fq.output,
        transcriptome="data/genome/extract_transcriptome/transcriptome.fa",
    output:
        res_dir + "result/align/transcriptome/run_minimap2_novel_transcriptome/{sample}.aligned.bam",
    params:
        align_map_bam_threads=config["align_map_bam_threads"],
    threads: config["align_map_bam_threads"]
    log:
        stdout="logs/align/transcriptome/run_minimap2_transcriptome/{sample}.stdout",
        stderr="logs/align/transcriptome/run_minimap2_transcriptome/{sample}.stderr",
    conda:
        "envs/minimap2.yaml"
    shell:
        """
        minimap2 -ax map-hifi -N 100 --sam-hit-only --for-only \
            -t {params.align_map_bam_threads}  \
            {input.transcriptome} {input.reads} > {output}
        """


rule add_BC_tags:
    input:
        res_dir + "result/align/transcriptome/run_minimap2_novel_transcriptome/{sample}.aligned.bam",
    output:
        res_dir + "result/align/transcriptome/add_BC_tags/{sample}.aligned.tagged.bam",
    log:
        stderr="logs/align/transcriptome/add_BC_tags/{sample}.stderr",
    conda:
        "envs/jvarkit.yaml"
    shell:
        """
        samjdk --samoutputformat BAM \
            -e '
            String s = record.getReadName(); 
            int u = s.indexOf(\\"_\\");record.setReadName(s.substring(u+1)); 
            record.setAttribute(\\"CB\\",s.substring(0,u));
            return record;' {input} > {output} 2> {log.stderr}
        """

rule merge_bam_tx:
    input:
        lambda wildcards: expand(
            res_dir + "result/align/transcriptome/add_BC_tags/{sample}.aligned.tagged.bam",
            sample=MERGE[wildcards.tissue]
        )
    output:
        res_dir + "result/align/transcriptome/merge/{tissue}.aligned.tagged.bam"
    log:
        stderr="logs/align/transcriptome/merge_bam1/{tissue}.stderr"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        if [ $(echo "{input}" | wc -w) -gt 1 ]; then
            samtools merge {output} {input} 2> {log.stderr}
            rm -rf {input}
        else
            mv {input} {output}
        fi
        """





##################### Align to genome #########################

rule align_minimap2_genome:
    input:
        reads=res_dir+"result/align/transcriptome/bam2fq/{sample}.fastq.gz",
        transcriptome="data/genome/gencode.v47.chr_patch_hapl_scaff.annotation.bed",
        genome="data/genome/GRCh38.p14.genome.fa",
    output:
        "/home/jiayiwang/loggedfs_david_penton/result/align/genome/run_minimap2_genome/{sample}.aligned.sorted.bam"
    params:
        align_map_bam_threads=config["align_map_bam_threads"],
        align_sort_bam_threads=config["align_sort_bam_threads"],
        align_sort_bam_memory_gb=config["align_sort_bam_memory_gb"],
    threads: config["align_sort_bam_threads"] + config["align_map_bam_threads"]
    log:
        stdout="logs/align/genome/run_minimap2_genome/{sample}.stdout",
        stderr="logs/align/genome/run_minimap2_genome/{sample}.stderr"
    conda:
        "envs/minimap2.yaml"
    shell:
        """
        minimap2 -ax splice:hq -uf --junc-bed {input.transcriptome} \
            -t {params.align_map_bam_threads}  \
            {input.genome} {input.reads} | samtools sort \
            -@ {params.align_map_bam_threads} \
            -m{params.align_sort_bam_memory_gb}g \
            -o {output} > {log.stdout} \
            2> {log.stderr}
        """

rule merge_bam_gm:
    input:
        lambda wildcards: expand(
            res_dir + "result/align/genome/run_minimap2_genome/{sample}.aligned.sorted.bam",
            sample=MERGE[wildcards.tissue]
        )
    output:
        res_dir + "result/align/genome/merge/{tissue}.aligned.tagged.bam"
    log:
        stderr="logs/align/genome/merge/{tissue}.stderr"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        if [ $(echo "{input}" | wc -w) -gt 1 ]; then
            samtools merge {output} {input} 2> {log.stderr}
            rm -rf {input}
        else
            mv {input} {output}
        fi
        """


rule add_BC_tags_genome:
    input:
        res_dir + "result/align/genome/merge/{tissue}.aligned.tagged.bam"
    output:
        "/home/jiayiwang/loggedfs_david_penton/result/align/genome/add_BC_tags/{tissue}.aligned.tagged.bam",
    log:
        stderr="logs/align/genome/add_BC_tags/{tissue}.stderr",
    conda:
        "envs/jvarkit.yaml"
    shell:
        """
        samjdk --samoutputformat BAM \
            -e '
            String s = record.getReadName(); 
            int u = s.indexOf(\\"_\\");record.setReadName(s.substring(u+1)); 
            record.setAttribute(\\"CB\\",s.substring(0,u));
            return record;' {input} > {output} 2> {log.stderr}
        """



######################## Quantification ##########################
rule sort_bams_by_barcode_transcriptome:
    input:
        res_dir+"result/align/transcriptome/merge/{tissue}.aligned.tagged.bam",
    output:
        res_dir+"result/quantify/transcriptome/sort_bams_by_barcode/{tissue}.aligned.sorted_barcode.bam",
    threads: config["sorting_threads"]
    params:
        sorting_threads=config["sorting_threads"],
        sorting_memory=config["sorting_memory"],
    log:
        "logs/quantify/transcriptome/sort_bams_by_barcode/{tissue}.out",
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools sort -@ {params.sorting_threads} -m{params.sorting_memory}g -t CB {input} -o {output} &> {log}
        """

# rule collate_bams:
#     input:
#         res_dir+"result/quantify/transcriptome/sort_bams_by_barcode/{tissue}.aligned.sorted_barcode.bam"
#     output:
#         res_dir+"result/quantify/transcriptome/collate_bams/{tissue}.aligned.sorted_barcode.collate.bam"
#     params:
#         sorting_memory=config["sorting_memory"]
#     threads: config["collating_threads"]
#     log:
#         "logs/quantify/transcriptome/collate_bams/{tissue}.log",
#     conda:
#         "envs/samtools.yaml"
#     shell:
#         """
#         samtools sort -n -@ {threads} -m{params.sorting_memory}g -o {output} {input}
#         """



rule oarfish_quantify:
    input:
        res_dir+"result/quantify/transcriptome/sort_bams_by_barcode/{tissue}.aligned.sorted_barcode.bam"
    output:
        res_dir+"result/quantify/transcriptome/oarfish/{tissue}.count.mtx",
        res_dir+"result/quantify/transcriptome/oarfish/{tissue}.barcodes.txt",
        res_dir+"result/quantify/transcriptome/oarfish/{tissue}.features.txt"
    log:
        "logs/quantify/transcriptome/oarfish_quantify/{tissue}.out"
    conda:
        "envs/oarfish.yaml"
    params:
        outdir=res_dir+"result/quantify/transcriptome/oarfish/{tissue}",
        thread=20
    shell:
        """
        oarfish --single-cell --model-coverage --filter-group no-filters --output {params.outdir} --alignments {input} --threads {params.thread} &> {log}
        """

rule quantify_extract_transcript_gene_mapping:
    input:
        "data/genome/gencode.v47.chr_patch_hapl_scaff.annotation.gtf",
    output:
        "data/extract_transcript_gene_mapping/transcript_gene_mapping.tsv",
    threads: 1
    log:
        "logs/quantify/transcriptome/extract_transcript_gene_mapping/log.out",
    conda:
        "envs/gffread.yaml"
    shell:
        """
        gffread {input} --table transcript_id,gene_id > {output} 2> {log};
        sed -i '1i \ttranscript_id\tgene_id' {output} 2>> {log}
        """

rule generate_sce:
    input:
        "scripts/r/generate_sce.R",
        tx="data/extract_transcript_gene_mapping/transcript_gene_mapping.tsv",
        mtx=res_dir+"result/quantify/transcriptome/oarfish/{tissue}.count.mtx",
        fts=res_dir+"result/quantify/transcriptome/oarfish/{tissue}.features.txt",
        bcd=res_dir+"result/quantify/transcriptome/oarfish/{tissue}.barcodes.txt"
    output:
        tsce=res_dir+"result/quantify/transcriptome/transcript_sce/{tissue}.rds",
        gsce=res_dir+"result/quantify/transcriptome/gene_sce/{tissue}.rds"
    log:
        "logs/quantify/transcriptome/generate_sce/{tissue}.out"
    shell:
        """
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        mtx={input.mtx} fts={input.fts} bcd={input.bcd} tx={input.tx} tsce={output.tsce} gsce={output.gsce}" {input[0]} {log}
        """

# rule quantify_split_bams_by_barcode_transcriptome:
#     input:
#          res_dir+"result/quantify/transcriptome/sort_bams_by_barcode/{tissue}.aligned.sorted_barcode.bam",
#     output:
#         directory(res_dir+"result/quantify/transcriptome/split_bams_by_barcode/{tissue}"),
#     params:
#         input_prefix="aligned.sorted_by_barcode",
#     threads: config["split_threads"]
#     log:
#         "logs/quantify/transcriptome/split_bams_by_barcode/{tissue}.out",
#     conda:
#         "envs/pysam.yaml"
#     script:
#         "scripts/py/split_barcodes.py"


# rule get_barcodes_transcriptome:
#     input:
#         rules.quantify_split_bams_by_barcode_transcriptome.output
#     output:
#         res_dir+"result/quantify/transcriptome/get_barcodes/{tissue}.barcodes.txt",
#     log:
#         "logs/quantify/transcriptome/get_barcodes/{tissue}.out",
#     conda:
#         "envs/pysam.yaml"
#     script:
#         "scripts/py/get_barcodes.py"


# rule quantify_collate_bams:
#     input:
#         indir= res_dir+"result/quantify/transcriptome/split_bams_by_barcode/{tissue}",
#         bc_list= res_dir+"result/quantify/transcriptome/get_barcodes/{tissue}.barcodes.txt",
#     output:
#         outdir=directory(res_dir+"result/quantify/transcriptome/collate_bams/{tissue}"),
#     params:
#         barcodes=lambda wc, input: [
#             line.strip() for line in open(input["bc_list"], "r") if line.strip()
#         ],
#         sorting_memory=config["sorting_memory"],
#     threads: config["collating_threads"]
#     log:
#         "logs/quantify/transcriptome/collate_bams/{tissue}.log",
#     conda:
#         "envs/samtools.yaml"
#     shell:
#         """
#         mkdir {output.outdir};
#         for barcode in {params.barcodes}; do
#             samtools sort -n -@ {threads} -m{params.sorting_memory}g -o {output}/aligned.collated.$barcode.bam {input.indir}/aligned.sorted_by_barcode_$barcode.bam
#         done
#         """


# rule quantify_run_oarfish:
#     input:
#         indir=res_dir+"result/quantify/transcriptome/collate_bams/{tissue}",
#         bc_list=res_dir+"result/quantify/transcriptome/get_barcodes/{tissue}.barcodes.txt",
#     output:
#         outdir=directory(res_dir+"result/quantify/transcriptome/run_oarfish/{tissue}"),
#     params:
#         barcodes=lambda wc, input: [
#             line.strip() for line in open(input["bc_list"], "r") if line.strip()
#         ],
#         n_bins=10,
#         filter_group="no-filters",
#     threads: config["quantify_threads"]
#     log:
#         "logs/quantify/transcriptome/run_oarfish/{tissue}.log",
#     conda:
#         "envs/oarfish.yaml"
#     shell:
#         """
#         for barcode in {params.barcodes}; do
#             oarfish --threads {threads} \
#                 --filter-group {params.filter_group} \
#                 --model-coverage \
#                 --alignments {input.indir}/aligned.collated.${{barcode}}.bam \
#                 --output {output.outdir}/${{barcode}}
#         done
#         """
# #--bins {params.n_bins} \



# rule quantify_build_transcriptome_sce:
#     input:
#         indir=res_dir+"result/quantify/transcriptome/run_oarfish/{tissue}",
#         bc_list=res_dir+"result/quantify/transcriptome/get_barcodes/{tissue}.barcodes.txt",
#     output:
#         output_path=res_dir+"result/quantify/transcriptome/quantify_build_transcriptome_sce/{tissue}/transcriptome_sce.rds",
#     params:
#         barcodes=lambda wc, input: [
#             line.strip() for line in open(input["bc_list"], "r") if line.strip()
#         ],
#         sample_id="{tissue}",
#     threads: config["quantify_threads"]
#     log:
#         "logs/quantify/transcriptome/quantify_build_transcriptome_sce/{tissue}/log.out",
#     conda:
#         "envs/merge_oarfish.yaml"
#     script:
#         "scripts/r/summarize_per_sample_transcriptome.R"





# rule quantify_build_genome_sce:
#     input:
#         indir=res_dir+"result/quantify/transcriptome/oarfish/{tissue}",
#         gene_transcriptome_mapping="data/genome/transcript_gene_mapping.tsv",
#         bc_list=res_dir+"result/quantify/transcriptome/get_barcodes/{tissue}.barcodes.txt",
#     output:
#         output_path=res_dir+"result/quantify/transcriptome/quantify_build_genome_sce/{tissue}/genome_sce.rds",
#     params:
#         barcodes=lambda wc, input: [
#             line.strip() for line in open(input["bc_list"], "r") if line.strip()
#         ],
#         sample_id="{tissue}",
#     threads: config["quantify_threads"]
#     log:
#         "logs/quantify/transcriptome/quantify_build_genome_sce/{tissue}/log.out",
#     conda:
#         "envs/merge_oarfish.yaml"
#     script:
#         "scripts/r/summarize_per_sample_genome.R"

################## Variant Calling #################

rule variant_calling_longshot:
    input:
        indir="results/variant/bam_by_barcodes/{merge}/",
        ref="data/genome/GRCh38.primary_assembly.genome.fa"
    output:
        outdir=directory("results/variant/longshot/{merge}")
    log:
        "logs/variant/longshot/{merge}.log.out"
    shell:
       """
        mkdir -p {output.outdir}
        for bam_file in {input.indir}/*.bam; do
            prefix=$(basename "$bam_file" .bam)
            samtools index "$bam_file"
            longshot --bam "$bam_file" --ref {input.ref} --out {output.outdir}/$prefix.vcf > {log} 2>&1
        done
       """

rule annotate_vcf:
    input:
        indir="results/variant/longshot/{merge}/"
    output:
        outdir=directory("results/variant/anno_vcf/{merge}")
    log:
        "logs/variant/anno_vcf/{merge}.log.out"
    shell:
        """
        mkdir -p {output.outdir}
        for vcf_file in {input.indir}/*.vcf; do
            prefix=$(basename "$vcf_file" .vcf)
            java -Xmx8g -jar ~/snpEff/snpEff.jar GRCh38.86 "$vcf_file" > {output.outdir}/$prefix.anno.vcf -s {output.outdir}/$prefix.html
        done
        """





# sinto filterbarcodes -b {sample}.aligned.sorted.bam -c barcodes.txt --outdir {sample}
# longshot --bam  GCCCAGTGACACGGTC.bam --ref ~/pacbio/data/genome/GRCh38.primary_assembly.genome.fa --out  GCCCAGTGACACGGTC.vcf
# java -Xmx8g -jar ~/snpEff/snpEff.jar GRCh38.86 GCCCAGTGACACGGTC.vcf > GCCCAGTGACACGGTC.anno.vcf


# to run it: snakemake --use-conda --use-singularity --cores 37 --singularity-args '-B /home/gibott/data,/home/gibott/pacbio_seq/exp_bulk/workflow/results,/home/gibott/illumina_seq/exp_sc/workflow/results' --rerun-incomplete
# mount the results folder of PacBio bulk for the extended bambu annotations
# mount the results folder of illumina single cell for the sce object
# use --rerun-incomplete if you have some completed samples in a step
