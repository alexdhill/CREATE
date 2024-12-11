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


process download_acc
{
    if (params.manage_resources)
    {
        cpus 1
        memory '1.GB'
    }
    errorStrategy "retry"
    maxRetries 10
    input:
        tuple(
            val("acc"),
            val("file"),
            val("md5")
        )
    output:
        tuple(
            val("${acc}"),
            path("*.fastq.gz")
        )
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Downloading reads..."
                echo "Sample: !{acc}"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            md5sum=""
            until [[ "!{md5}" == "${md5sum}" ]]; do
                until `wget !{file}`; do
                    echo "Download failed... retrying in 5s"
                    sleep 5
                done
                md5sum="$(md5sum $(basename !{file}) | awk '{print $1}')"
            done
        '''
}