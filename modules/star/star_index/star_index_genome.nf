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
 

process star_index_genome
{
    publishDir "${params.outdir}/", mode: 'copy', overwrite: params.force
    if (params.manage_resources)
    {
        cpus 32
        memory '256.GB' // TODO
    }
    input:
        path(genome)
    output:
        path("${params.genome}v${params.genome=='T2T'?'2':params.version}_discover_index_v*.star")
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Creating complete STAR index"
                echo "Genome: !{genome}"
            fi
            version="!{params.version}"
            if [[ "!{params.genome}" == "T2T" ]]; then
                version="2"
            fi
            if [[ "!{params.log}" == "DEBUG" ]]; then
                set -x
            fi

            zcat !{genome} > genome.fa
            STAR --runMode genomeGenerate \
                --genomeDir !{params.genome}v${version}_discover_index_v$(STAR --version).star \
                --genomeFastaFiles genome.fa \
                --runThreadN 32
        '''
}
