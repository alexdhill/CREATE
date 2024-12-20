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


process gather_ftp
{
    if (params.manage_resources)
    {
        cpus 1
        memory '1.GB'
    }
    executor = "local"
    input:
        path(acc)
    output:
        path("samples.url")
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Collecting FTP links..."
                echo "Sample: !{acc}"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            for sra in `xargs < !{acc}`; do
                until `ffq --ftp ${sra} >> samples.json`; do
                    echo "Failed at ${sra}, retrying in 5s..."
                    sleep 5
                done
            done

            jq '.[] | "\\(.accession),\\(.url),\\(.md5)"' samples.json \
            | sed 's/"//g' \
            > samples.url
        '''
}