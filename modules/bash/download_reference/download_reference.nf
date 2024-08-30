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
 

process download_reference
{
    publishDir "${params.outdir}", mode: 'copy', enable: params.keep, overwrite: params.force
    if (params.manage_resources)
    {
        cpus 1
        memory '1.GB'
    }
    output:
        path("*_genome.fa.gz")
    shell:
        if (params.genome=="T2T")
        {
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Downloading T2Tv2 genome..."
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            wget -qO- https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz \
            | gzip -cd \
            | sed 's/\\s.*//' \
            | gzip --fast \
            > T2Tv2_genome.fa.gz
        '''
        }
        else
        {
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Downloading HG38 genome..."
            fi
            if [[ !{params.log}=="DEBUG" ]]; then
                set -x
            fi

            wget -qO- ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_!{params.version}/GRCh38.primary_assembly.genome.fa.gz \
            | gzip -cd \
            | sed 's/\\s.*//' \
            | gzip --fast \
            > HG38v!{params.version}_genome.fa.gz
        '''
        }
}