data = "data/fastq/{id}.fq"
ids, = glob_wildcards(data)

configfile: "config.yaml"

# rule all:
#     input: 
#         expand("result/{id}_subsample_{subsample}_Aligned.out.sam", id=ids, subsample=config["subsample"])

rule get_subsample:
    input:
        data
    output:
        temp("{id}_tmp_{subsample}.fq")
    params:
        "{subsample}"
    shell:
        "head -n {params} {input} > {output}"

rule mapping:
    input:
        fastq="{id}_tmp_{subsample}.fq",
        index=lambda wc: config["sample_dict"][wc.id]
    output: 
        "result/{id}_subsample_{subsample}_Aligned.out.sam"
    params:
        "{subsample}"
    shell:
        "STAR "
        "--genomeDir {input.index} "
        "--readFilesIn {input.fastq} "
        "--outFileNamePrefix result/{wildcards.id}_subsample_{params}_"