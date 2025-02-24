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
 * connor was here 
 */
 

process trim_reads_np
{
    publishDir "${params.outdir}/reads/trimmed/", mode: 'copy', enable: params.keep, overwrite: params.force
    if (params.manage_resources)
    {
        cpus 8
        memory '16.GB'
    }
    input:
        tuple(
            val(sample),
            val(nreads),
            path(read)
        )
    output:
        tuple(
            val("${sample}"),
            env(NREADS),
            path("${sample}_filtered_trimmed.fq.gz")
        )
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Trimming nanopore reads..."
                echo "Sample: !{sample}"
                echo "Read: !{read}"
                echo "n Reads: !{nreads}"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            chopper \
                -i !{read} \
                -q 12 \
                -l 70 \
                --threads !{task.cpus} \
            | pigz \
            > !{sample}_filtered_trimmed.fq.gz

            NREADS=`gzip -cd !{sample}_filtered_trimmed.fq.gz \
            | wc -l \
            | awk '{print $1/4}'`
        '''
}
