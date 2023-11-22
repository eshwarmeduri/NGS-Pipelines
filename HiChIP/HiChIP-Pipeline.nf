#!/usr/bin/env nextflow -DSL2

params.reads = "$baseDir/Files/*.r_{1,2}.fq.gz"
genome="/tmp/hg38.fa"
chrSizes="/tmp/hg38.chrom.sizes"
params.cores=24


workflow {
    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true )
    
    trimmed_ch=trimGalore(read_pairs_ch)

    genome_ch = mapgenome(trimmed_ch)

    ligat_ch = validligation(genome_ch)

    remdup_ch=remDupGenBAM(ligat_ch)


}

process trimGalore {
    tag "cutadapt on $sample_id"

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

    tag "BWA alignment sample ${sample_id}"

    publishDir "Output", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*")

    script:
    """
    bwa mem -5SP -T0 -t $params.cores $genome ${sample_id}.r_1_val_1.fq.gz ${sample_id}.r_2_val_2.fq.gz -o ${sample_id}.sam
    """

}

process validligation {

    tag "Recording valid ligation events for sample ${sample_id}"

    publishDir "Output", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*pairsam")

    script:
    """
    pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in $params.cores --nproc-out $params.cores --chroms-path $chrSizes ${sample_id}.sam >  ${sample_id}.pairsam
    pairtools sort --nproc $params.cores --tmpdir=$baseDir/ ${sample_id}.pairsam > ${sample_id}.sorted.pairsam
    
    """

}

process remDupGenBAM {
    tag "Removing duplicates and generating BAM files for ${sample_id}"

    publishDir "Output", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(reads_file)

    output:
    tuple val(sample_id), path("*pairsam")
    tuple val(sample_id), path("*dedup.stats.txt")
    tuple val(sample_id), path("*bam")

    script:
    """
    pairtools dedup --nproc-in $params.cores --nproc-out $params.cores --mark-dups --output-stats ${sample_id}.dedup.stats.txt --output ${sample_id}.dedup.pairsam ${sample_id}.sorted.pairsam
    pairtools split --nproc-in $params.cores --nproc-out $params.cores --output-pairs ${sample_id}.mapped.pairs --output-sam ${sample_id}.unsorted.bam ${sample_id}.dedup.pairsam
    samtools sort -@ $params.cores -o ${sample_id}.mapped.PT.bam ${sample_id}.unsorted.bam

    """
}
