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

process flair_collapse
{
    publishDir "${params.outdir}/isoforms", mode: 'copy', overwrite: params.force, enable: params.keep
    if (params.manage_resources)
    {
        cpus 8
        memory '64.GB' // TODO
    }
    input:
    tuple(
            path(region),
            path(reads),
            path(reference)
    )
    output:
        tuple(
            path("*.isoforms.fa"),
            path("*.isoforms.bed"),
            path("*.isoforms.gtf"),
            path("*.isoform.read.map.txt")
        )
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Running FLAIR discovery"
                echo "BEDs:\n!{region}"
                echo "Reads:\n!{reads}"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            base=$(basename -s .bed !{region})

            gzip -cd !{reference}/*genome.fa.gz > genome.fa
            gzip -cd !{reference}/*_complete_annotation.gtf.gz > annotation.gtf
            flair collapse \
                --genome genome.fa \
                --query !{region} \
                --reads !{reads} \
                --gtf annotation.gtf \
                --annotation_reliant generate \
                --stringent \
                --check_splice \
                --generate_map \
                --threads !{task.cpus} \
                --output $base
        '''
}