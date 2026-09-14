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


process multiqc_report_short
{
    publishDir "${params.outdir}/report/", mode: 'copy', overwrite: params.force
    container 'alexdhill/create:multiqc-1.25'
    conda projectDir+'/bin/conda/modules/multiqc.yaml'
    if (params.manage_resources)
    {
        cpus 2
        memory '16.GB'
        time '2h'
    }
    input:
        tuple(
            path(read_reports),
            path(quants),
            path(parameters)
        )
    output:
        tuple(
            path("multiqc_report.html"),
            path("multiqc_report_data/")
        )
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Summarizing short read quantification"
                echo "Read reports: $(ls -d !{read_reports} | wc -l)"
                echo "Quantifications: $(ls -d !{quants} | wc -l)"
                echo "User parameters: !{parameters}"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            params="--force --interactive"
            if [[ "!{parameters}" != "NULL" ]]; then
                params="$(jq '.multiqc | to_entries | .[] | "\\(.key)=\\(.value)"' !{parameters} | xargs | sed 's/=true//g')"
            fi

            multiqc ${params} \
                --filename multiqc_report.html \
                --outdir . \
                .

            if [[ ! -e multiqc_report.html ]]; then
                echo "\033[1;31mERR: MultiQC report failed\033[0m" 1>&2
                exit 1
            fi
        '''
}
