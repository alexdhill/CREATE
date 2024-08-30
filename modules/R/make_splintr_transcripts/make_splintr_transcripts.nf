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


process make_splintr_transcripts
{
    publishDir "${params.outdir}/", mode: 'copy', enable: params.keep, overwrite: params.force
    if (params.manage_resources)
    {
        cpus 1
        memory '24.GB'
    }
    input:
        tuple(
            path(annotation),
            path(genome)
        )
    output:
        tuple(
            path("*_splintr_transcripts.fa.gz"),
            path("*_splintr_map.tx2g"),
            path("*_splintr_map.tx3g")
        )
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Generating spliced+intron transcripts"
                echo "Annotation: !{annotation}"
                echo "Genome: !{genome}"
            fi
            version="!{params.version}"
            if [[ "!{params.genome}"=="T2T" ]]; then
                version="2"
            fi
            verbose=""
            if [[ "!{params.log}" == "DEBUG" ]]; then
                verbose="--verbose"
                set -x
            fi

            zcat !{genome} > genome.fa
            zcat !{annotation} > genes.gtf
            Rscript ${verbose} !{projectDir}/bin/R/make_splintr_transcripts.R \
                -a genes.gtf -g genome.fa \
                -v !{params.genome}v${version}
        '''
}