/*
 * REQUIRED NOTICE: Copyright (c) 2020-2023, Regents of the University of California
 * All rights reserved. https://polyformproject.org/licenses/noncommercial/1.0.0
 * 
 * This software was developed by the Daniel Kim lab at the University of California, Santa Cruz.
 * Authors: Roman E. Reggiardo, Vikas Peddu, Alex D. Hill
 * 
 * The licensor grants you a copyright license for the software to do everything you might do with
 * the software that would otherwise infringe the licensor’s copyright in it for any permitted
 * purpose.
 * 
 * As far as the law allows, the software comes as is, without any warranty or condition, and the
 * licensor will not be liable to you for any damages arising out of these terms or the use or
 * nature of the software, under any kind of legal claim.
 */
 

process star_align_paired
{
    publishDir "${params.outdir}/align", mode: 'copy', enabled: params.keep, overwrite: params.force
    if (params.manage_resources)
    {
        cpus 8
        memory '64.GB' // TODO
    }
    input:
        tuple(
            val(sample),
            path(read_1),
            path(read_2),
            val(nreads),
            path(reference)
        )
    output:
        tuple(
            path("${sample}.bam"),
            path("${sample}_unmapped.bam")
        )
    shell:
        '''
            if [[ !{params.log}=="INFO" || !{params.log}=="DEBUG" ]]; then
                echo "Creating complete salmon index"
                echo "Sample: !{sample}"
                echo "Read 1: !{read_1}"
                echo "Read 2: !{read_2}"
                echo "Reference: !{reference}"
            fi
            if [[ !{params.log}=="DEBUG" ]]; then
                set -x
            fi

            STAR \
                --genomeDir !{reference}/*.star \
                --readFilesIn !{read_1} !{read_2} \
                --quantMode TranscriptomeSAM \
                --outFileNamePrefix !{sample}. \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMunmapped Within \
                --runThreadN 8

            if [[ ! -e !{sample}.Aligned.sortedByCoord.out.bam ]]; then
                echo "\033[[1;31mERR: STAR alignment failed\033[0m" 1>&2
                exit 1
            fi

            samtools view -f 4 !{sample}.Aligned.sortedByCoord.out.bam > !{sample}_unmapped.bam &
            samtools view -F 4 !{sample}.Aligned.sortedByCoord.out.bam > !{sample}.bam
        '''
}