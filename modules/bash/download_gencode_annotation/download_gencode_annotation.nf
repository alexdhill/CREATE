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


process download_gencode_annotation
{
    publishDir "${params.outdir}/", mode: 'copy', enabled: params.keep, overwrite: params.force
    if (params.manage_resources)
    {
        cpus 1
        memory '1.GB'
    }
    output:
        path("${params.genome}v${params.genome=='T2T'?'2':params.version}_gencode_annotation.gtf.gz")
    shell:
        if (params.isoquant)
        {
            '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Preparing isquant annotation..."
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                echo "Isoquant: !{params.isoquant}"
            fi
            '''
        }
        else if (params.genome=="T2T")
        {
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Downloading T2T annotation..."
                echo "... using T2Tv2"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            wget -qO- 'https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13.draft_v2.0.gene_annotation.gff3' \
            | sed \
                -e's/transcript_id=/new_transcript_id=/' \
                -e's/source_transcript=/transcript_id=/' \
                -e's/gene_id=/new_gene_id=/' \
                -e's/source_gene=/gene_id=/' \
                -e's/;/; /g' \
                -e's/=/ /g' \
            | grep -v 'StringTie' \
            | grep -v 'gene_biotype unknown' \
            | grep -v "^#" \
            | gzip --best \
            > T2Tv2_gencode_annotation.gtf.gz
        '''
        }
        else
        {
        '''
            #!/bin/bash

            if [[ !{params.log}=="INFO" || !{params.log}=="DEBUG" ]]; then
                echo "Downloading HG38 annotation"
                echo "... using HG38v!{params.version}"
            fi
            if [[ !{params.log}=="DEBUG" ]]; then
                set -x
            fi

            wget -qO- ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_!{params.version}/gencode.v!{params.version}.annotation.gtf.gz \
            | gzip -cd \
            | sed 's/|/ /g' \
            | gzip --best \
            > HG38v!{params.version}_gencode_annotation.gtf.gz
        '''
        }
}