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
 

process star_index
{
    publishDir "${params.outdir}/", mode: 'copy', overwrite: params.force
    if (params.manage_resources)
    {
        cpus 8
        memory '64.GB' // TODO
    }
    input:
        tuple(
            path(annotation),
            path(genome)
        )
    output:
        path("*.star")
    shell:
        '''
            if [[ !{params.log}=="INFO" || !{params.log}=="DEBUG" ]]; then
                echo "Creating complete STAR index"
                echo "Genome: !{genome}"
                echo "Annotation: !{annotation}"
            fi
            version="!{params.version}"
            if [[ "!{params.genome}"=="T2T" ]]; then
                version="2"
            fi
            if [[ !{params.log}=="DEBUG" ]]; then
                set -x
            fi

            zcat !{annotation} > annotation.gtf
            zcat !{genome} > genome.fa

            STAR --runMode genomeGenerate \
                --genomeDir !{params.genome}v${version}_index_v$(STAR --version).star \
                --genomeFastaFiles genome.fa \
                --sjdbGTFfile annotation.gtf \
                --sjdbOverhang 149 \
                --limitSjdbInsertNsj 6000000 \
                --runThreadN !{task.cpus}
        '''
}