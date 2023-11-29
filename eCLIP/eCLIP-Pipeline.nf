#!/usr/bin/env nextflow -DSL2

params.reads = "$baseDir/Files/*.r_{1,2}.fq.gz"
repbase="/tmp/STAR-GENOMES-hg38-RepBase"
transcriptome="/tmp/STAR-GENOMES-GRCh38.gencode-RL50"
params.cores=1

workflow {
    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true )
    
    trimmed_ch=cutadapt(read_pairs_ch)

    genome_ch = stargenome(trimmed_ch)

    remdup_ch = remdup(genome_ch)


}

process cutadapt {
    tag "cutadapt on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path ("*_{1,2}.TrTr.fq.gz")


    script:
    """

    cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6 -m 18 -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CTTCCGATCTACAAGTT  -g CTTCCGATCTTGGTCCT -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA  -A CTTGTAGATCGGAAG  -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG  -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT  -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT -o ${sample_id}_1.Tr.fq.gz -p ${sample_id}_2.Tr.fq.gz ${reads[0]} ${reads[1]}

    cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 5 --quality-cutoff 6 -m 18 -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT -o ${sample_id}_1.TrTr.fq.gz -p ${sample_id}_2.TrTr.fq.gz ${sample_id}_1.Tr.fq.gz ${sample_id}_2.Tr.fq.gz

    """

}

process stargenome {

    tag "STAR align for ${sample_id}_1.TrTr.fq.gz ${sample_id}_2.TrTr.fq.gz against ${repbase} "

    publishDir "Output/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path ("*")

    script:
    """
    STAR --runMode alignReads --runThreadN $params.cores --genomeDir ${repbase} --genomeLoad NoSharedMemory --alignEndsType EndToEnd --outSAMunmapped Within --outFilterMultimapNmax 30 --outFilterMultimapScoreRange 1 --outFileNamePrefix ${sample_id}.sorted.repbase. --outSAMtype BAM Unsorted --outFilterType BySJout --outBAMcompression 10 --outReadsUnmapped Fastx --outFilterScoreMin 10 --outSAMattrRGline ID:foo --outSAMattributes All --outSAMmode Full --outStd Log --readFilesCommand zcat --readFilesIn ${sample_id}_1.TrTr.fq.gz ${sample_id}_2.TrTr.fq.gz

    mv ${sample_id}.sorted.repbase.*1 ${sample_id}.repbase.unmapped.r_1.fq
    mv ${sample_id}.sorted.repbase.*2 ${sample_id}.repbase.unmapped.r_2.fq

    STAR --runMode alignReads --runThreadN $params.cores --genomeDir ${transcriptome} --genomeLoad NoSharedMemory --readFilesIn ${sample_id}.repbase.unmapped.r_1.fq ${sample_id}.repbase.unmapped.r_2.fq --outFileNamePrefix ${sample_id}.genome.sorted.STAR. --outSAMattributes All --outSAMtype BAM SortedByCoordinate --alignEndsType EndToEnd --outFilterMultimapNmax 1 --outFilterScoreMin 10

    """

}

process remdup {

    tag "Remove PCR duplicates for sample ${sample_id}"

    publishDir "Output/${sample_id}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*")

    script:
    """
    java -jar /usr/picard/picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT INPUT=${sample_id}.genome.sorted.STAR.Aligned.sortedByCoord.out.bam OUTPUT=${sample_id}.genome.STAR.dedup.bam METRICS_FILE=${sample_id}.dedup_metrics.txt REMOVE_DUPLICATES=true CREATE_INDEX=true

    """

}
