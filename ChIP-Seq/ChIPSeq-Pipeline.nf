#!/usr/bin/env nextflow -DSL2

params.reads = "$baseDir/Files/*.r_{1,2}.fq.gz" //$baseDir is the current working directory
refGenome="/tmp/hg38" //hg38 reference genome
params.cores=1

workflow {
    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true )
    
    trimmed_ch=trimGalore(read_pairs_ch) //trim adapters and FastQC

    genome_ch = mapgenome(trimmed_ch) // mapping against hg38 genome using Bowtie2 and sortBAM

    remdup_ch = remdup(genome_ch)  // Removes PCR duplicates

    convertBW_ch = convertBigWig(remdup_ch) // converts BAM to BigWig

}

process trimGalore {
    tag "cutadapt on $sample_id"

    publishDir "Output/${sample_id}", mode: 'copy', pattern: "*.html", overwrite: true

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path ("*.fq.gz")

    script:
    """
    trim_galore --fastqc --cores $params.cores --paired ${reads[0]} ${reads[1]}

    """
}


process mapgenome {

    tag "Bowtie2 alignment sample ${sample_id}"

    publishDir "Output/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*")

    script:
    """
    bowtie2 -p $params.cores -x ${refGenome} -1 ${sample_id}.r_1_val_1.fq.gz -2 ${sample_id}.r_2_val_2.fq.gz | grep -v "XS:i:" | awk '{if(\$2==147 || \$2==83 || \$2==99 || \$2==163 || \$2==81 || \$2==97 || \$2==145 || \$2==161 || \$1~/^@/) print \$0}' | samtools view -bSF4 - > ${sample_id}.bam
    samtools sort ${sample_id}.bam -o ${sample_id}.sorted.bam
    rm ${sample_id}.bam

    """
}

process remdup {

    tag "Remove PCR duplicates for sample ${sample_id}"

    publishDir "Output/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.dedup.*")

    script:
    """
    java -jar /usr/picard/picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT INPUT=${sample_id}.sorted.bam OUTPUT=${sample_id}.dedup.bam METRICS_FILE=${sample_id}.dedup_metrics.txt REMOVE_DUPLICATES=true CREATE_INDEX=true

    """
}

process convertBigWig {

    tag "Convert to BigWig for ${sample_id}"

    publishDir "Output/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.dedup.*")

    script:
    """
    bamCoverage --bam ${sample_id}.dedup.bam -o ${sample_id}.dedup.CPM.bw --binSize 20 --numberOfProcessors $params.cores -of bigwig --extendReads --normalizeUsing CPM

    """

}
