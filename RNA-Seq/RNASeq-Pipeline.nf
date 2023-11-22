#!/usr/bin/env nextflow -DSL2

params.reads = "$baseDir/Files/*.r_{1,2}.fq.gz" //$baseDir is the current working directory
transcriptome="/tmp/STAR-GENOMES-GRCh38.gencode-RL50" //hg38 reference genome for read length 50. 
params.cores=1
params.gtfFile='/tmp/gencode.v26.primary_assembly.annotation.gtf' //GTF file for feature counts

workflow {
    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true )
    
    trimmed_ch=trimGalore(read_pairs_ch) //trim adapters and FastQC

    genome_ch = stargenome(trimmed_ch) // mapping against hg38 genome using STAR

    remdup_ch = remdup(genome_ch) // Removes PCR duplicates

    countreads_ch=countReads(remdup_ch) // Get Featurecounts for mapped BAM files

    convertBW_ch = convertBigWig(remdup_ch) // converts BAM to BigWig


}

process trimGalore {
    tag "cutadapt on $sample_id"

    publishDir "Output/${sample_id}", mode: 'copy', pattern: "*.html", overwrite: true

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path ("*")


    script:
    """
    trim_galore --fastqc --cores $params.cores --paired ${reads[0]} ${reads[1]}

    """

}

process stargenome {

    tag "STAR alignment sample ${sample_id}"

    publishDir "Output/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.STAR.*")

    script:
    """

    STAR --genomeDir ${transcriptome} --outFileNamePrefix ${sample_id}.STAR. --readFilesIn ${sample_id}.r_1_val_1.fq.gz ${sample_id}.r_2_val_2.fq.gz --alignEndsType EndToEnd --runThreadN $params.cores --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat

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
    java -jar /usr/picard/picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT INPUT=${sample_id}.STAR.Aligned.sortedByCoord.out.bam OUTPUT=${sample_id}.STAR.dedup.bam METRICS_FILE=${sample_id}.dedup_metrics.txt REMOVE_DUPLICATES=true CREATE_INDEX=true

    """

}

process countReads {
    tag "Getting featurecounts for ${sample_id}"

    publishDir "Output/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path ("*.featureCounts.*")

    script:
    """
    featureCounts -F GTF -a $params.gtfFile -B -p -o ${sample_id}.featureCounts.txt --primary --ignoreDup ${sample_id}.STAR.dedup.bam

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
    bamCoverage -b ${sample_id}.STAR.dedup.bam --numberOfProcessors 24 --effectiveGenomeSize 2701495761 --normalizeUsing CPM -o ${sample_id}.dedup.CPM.bw
    """

}
