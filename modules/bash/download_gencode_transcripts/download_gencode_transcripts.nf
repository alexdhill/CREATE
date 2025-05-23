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
 

process download_gencode_transcripts
{
    publishDir "${params.outdir}/", mode: 'copy', enabled: params.keep, overwrite: params.force
    if (params.manage_resources)
    {
        cpus 4
        memory '4.GB'
    }
    input:
        tuple(
            path(annotation),
            path(genome)
        )
    output:
        path("${params.genome}v${params.genome=='T2T'?'2':params.version}_gencode_transcripts.fa.gz")
    shell:
        if (params.genome=="T2T")
        {
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Generating T2Tv2 transcripts..."
                echo "Annotation: !{annotation}"
                echo "Genome: !{genome}"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            pigz -cdp !{task.cpus} !{annotation} > genes
            pigz -cdp !{task.cpus} !{genome} > genome
            gffread -w txome -g genome genes
            pigz -cp !{task.cpus} txome \
            > T2Tv2_gencode_transcripts.fa.gz
        '''
        }
        else
        {
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Downloading HG38v!{params.version} transcripts..."
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            wget -qO- ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_!{params.version}/gencode.v!{params.version}.transcripts.fa.gz \
            > HG38v!{params.version}_gencode_transcripts.fa.gz
        '''
        }
}