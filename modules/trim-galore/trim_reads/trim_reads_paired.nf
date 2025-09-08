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
 

process trim_reads_paired
{
    publishDir "${params.outdir}/reads/trimmed", mode: 'copy', enabled: params.keep, overwrite: params.force
    if (params.manage_resources)
    {
        cpus 8
        memory '16.GB'
        time '30m'
    }
    input:
        tuple(
            val(sample),
            val(nreads),
            path(read_1),
            path(read_2),
            path(parameters)
        )
    output:
        tuple(
            val("${sample}"),
            path("${sample}_val_1.fq.gz"),
            path("${sample}_val_2.fq.gz"),
            env(NREADS)
        )
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Trimming paired reads..."
                echo "Sample: !{sample}"
                echo "Read 1: !{read_1}"
                echo "Read 2: !{read_2}"
                echo "n Reads: !{nreads}"
                echo "User parameters: !{parameters}"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi
            
            params="-l 75 --2color 20"
            if [[ "${parameters}" != "NULL" ]]; then
                params="$(jq '."trim-galore" | to_entries | .[] | "\\(.key)=\\(.value)"' flags.json | xargs | sed 's/=true//g')"
            fi

            trim_galore --paired --gzip  !{read_1} !{read_2} \
                --basename !{sample} \
                -j !{task.cpus} --output_dir . \
                ${params}

            NREADS=`gzip -cd !{sample}_val_1.fq.gz \
            | awk 'END {print NR/4}'`
        '''
}
