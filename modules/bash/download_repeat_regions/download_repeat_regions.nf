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
 

process download_repeat_regions
{
    publishDir "${params.outdir}/", mode: 'copy', enabled: params.keep, overwrite: params.force
    if (params.manage_resources)
    {
        cpus 1
        memory '1.GB'
    }
    output:
        path("*_repeat_regions.bed")
    shell:
        if (params.genome=="T2T")
        {
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Downloading T2T repeat regions..."
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            wget -qO- https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed \
            > T2Tv2_repeat_regions.bed
        '''
        }
        else
        {
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Downloading !{params.genome} repeat regions..."
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            build="hg38"
            if [[ "$(echo !{} | sed 's/\\B\\w*//' )" == "M" ]]; then
                build="mm39"
            fi

            bash !{projectDir}/bin/sql/repeatmasker_regions.mysql ${build} > raw_rmsk.bed
            for i in `seq 1 22` X Y; do
                grep -w "^chr$i" raw_rmsk.bed >> !{params.genome}v!{params.version}_repeat_regions.bed
            done
        '''
        }
}