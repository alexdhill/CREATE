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


process fastqc_report_paired
{
    publishDir "${params.outdir}/report/fastqc/${stage}/", mode: 'copy', enabled: params.keep, overwrite: params.force
    container 'alexdhill/create:fastqc-0.12.1'
    conda projectDir+'/bin/conda/modules/fastqc.yaml'
    if (params.manage_resources)
    {
        cpus 4
        memory '8.GB'
        time '2h'
    }
    input:
        tuple(
            val(sample),
            val(stage),
            path(read_1),
            path(read_2),
            path(parameters)
        )
    output:
        path("${sample}_${stage}/")
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Reporting paired read quality"
                echo "Sample: !{sample} (!{stage})"
                echo "Read 1: !{read_1}"
                echo "Read 2: !{read_2}"
                echo "User parameters: !{parameters}"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            params="--quiet --nogroup"
            if [[ "!{parameters}" != "NULL" ]]; then
                params="$(jq '.fastqc | to_entries | .[] | "\\(.key)=\\(.value)"' !{parameters} | xargs | sed 's/=true//g')"
            fi

            mkdir -p !{sample}_!{stage}

            fastqc ${params} \
                -t !{task.cpus} \
                --outdir !{sample}_!{stage} \
                !{read_1} !{read_2}

            if [[ $(ls !{sample}_!{stage}/*_fastqc.zip | wc -l) -ne 2 ]]; then
                echo "\033[1;31mERR: FastQC report failed\033[0m" 1>&2
                exit 1
            fi
        '''
}
