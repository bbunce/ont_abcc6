SAMPLES = ["07", "08", "09", "10", "11"]
#READS = [x for x in range(27)]
READS = [0,1,2]

rule align:
    input:
        "../data/ref/chr16.fa.gz",
        "../data/samples/ADC563_pass_barcode{sample}_522c8849_{read}.fastq.gz",
    output:
        "../aligned/sam/{sample}_{read}.sam"
    shell:
        "minimap2 -ax map-ont {input} > {output}"

rule sam_to_bam:
    input:
        "../aligned/sam/{sample}_{read}.sam"
    output:
        "../aligned/bam/{sample}_{read}.bam"
    shell:
        "samtools view -S -b {input} > {output}"

#todo this rule is not correct - it merges all files into the one
rule merge_bam:
    input:
        expand("../aligned/bam/{sample}_{read}.bam", sample=SAMPLES, read=READS)
    output:
        "../aligned/bam/merged/{sample}_merged.bam"
    shell:
        "samtools merge {output} {input}"

rule bam_sort:
    input:
        "../aligned/bam/merged/{sample}_merged.bam"
    output:
        "../aligned/sorted/{sample}_merged.sorted.bam"
    shell:
        "samtools sort {input} -o {output}"

rule bam_index:
    input:
        "../aligned/sorted/{sample}_merged.sorted.bam"
    output:
        "../aligned/sorted/{sample}_merged.sorted.bam.bai"
    shell:
        "samtools index {input}"