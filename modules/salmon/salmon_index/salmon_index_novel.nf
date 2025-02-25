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
 

process salmon_index_novel
{
    publishDir "${params.dump}/", mode: 'copy', overwrite: params.force, enable: params.dump!=''
    if (params.manage_resources)
    {
        cpus 8
        memory '24.GB'
    }
    input:
        path(transcripts)
    output:
        path("*.sidx")
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Creating complete salmon index"
                echo "Transcripts: !{transcripts}"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            pigz -cdp !{task.cpus} !{transcripts} \
            > transcripts.fa

            salmon index -t transcripts.fa \
                -p !{task.cpus} \
                -i novel_short_index_v$(salmon --version | awk '{print $2}').sidx
        '''
}