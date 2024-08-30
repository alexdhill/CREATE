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
 

process alevin_align
{
    publishDir "${params.outdir}/", mode: 'copy', enable: params.keep, overwrite: params.force
    if (params.manage_resources) {
        cpus 8
    }
    input:
        tuple(
            val(sample),
            path(read1),
            path(read2),
            path(reference)
        )
    output:
        tuple(
            val("${sample}"),
            path("alevin/${sample}")
        )
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Aligning reads..."
                echo "Sample: !{sample}"
                echo "Read 1: !{read1}"
                echo "Read 2: !{read2}"
                echo "Reference: !{reference}"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            mkdir -p alevin
            salmon alevin -l ISR \
                --!{params.chemistry} \
                -1 !{read1} -2 !{read2} \
                -i !{reference}/*splintr_index* \
                --sketch -p !{task.cpus} \
                -o alevin/!{sample}
        '''
}