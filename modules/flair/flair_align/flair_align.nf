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
 

process flair_align
{
    publishDir "${params.outdir}/align/nanopore", mode: 'copy', overwrite: params.force, enable: params.keep
    if (params.manage_resources)
    {
        cpus 8
        memory '64.GB' // TODO
    }
    input:
        tuple(
            val(sample),
            val(nreads),
            path(read),
            path(reference)
        )
    output:
        tuple(
            val("${sample}"),
            val("${nreads}"),
            path("${sample}.bed")
        )
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Running FLAIR alignment"
                echo "Sample: !{sample} (!{nreads} reads)"
                echo "Read: !{read}"
                echo "Reference: !{reference}"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            flair align \
                --genome !{reference}/*_genome.fa.gz \
                --reads !{read} \
                --threads 8 \
                --output !{sample}
        '''
}