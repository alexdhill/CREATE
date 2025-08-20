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
 

process minimap2_index
{
    publishDir "${params.outdir}/", mode: 'copy', overwrite: params.force
    if (params.manage_resources)
    {
        cpus 3
        memory '16.GB'
    }
    input:
        path(transcripts)
    output:
        path("*.mmi")
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Creating complete Minimap2 Index"
                echo "Transcripts: !{transcripts}"
            fi
            version="!{params.version}"
            if [[ "!{params.genome}" == "T2T" ]]; then
                version="2"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            mkfifo transcripts
            pigz -cdp !{task.cpus} !{transcripts} > transcripts &
            minimap2 transcripts \
                -t 3 \
                -d !{params.genome}v${version}_long_index_v$(minimap2 --version | awk -F'-' '{print $1}').mmi
        '''
}