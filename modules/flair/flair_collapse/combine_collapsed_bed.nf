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

process combine_collapsed_bed
{
    publishDir "${params.outdir}/combine_collapsed_bed", mode: 'copy', overwrite: params.force, enable: params.keep
    if (params.manage_resources)
    {
        cpus 1
        memory '16.GB' // TODO
    }
    //for some reason this only works with file and not path
    input:
            tuple val(fasta), file(isoform_fastas)
            tuple val(bed), file(isoform_beds)
            tuple val(gtf), file(isoform_gtfs)
            tuple val(readmap), file(isoform_read_maps)
    output:
        tuple(
            path("collapse_combined.fasta"),
            path("collapse_combined.bed"),
            path("collapse_combined.gtf"),
            path("collapse_combined.read.map.txt")
        )
        /*        tuple(
            tuple(val('fasta'), path("collapse_combined.fasta")),
            tuple(val('bed'),path("collapse_combined.bed")),
            tuple(val('gtf'),path("collapse_combined.gtf")),
            tuple(val('readmap'),path("collapse_combined.read.map.txt"))
        )*/
    shell:
        '''
            if [[ "!{params.log}" == "INFO" || "!{params.log}" == "DEBUG" ]]; then
                echo "Combining per-chromosome FLAIR collapsed bed files"
                echo "fastas:\n!{isoform_fastas}"
                echo "beds:\n!{isoform_beds}"
                echo "gtfs:\n!{isoform_gtfs}"
                echo "readmaps:\n!{isoform_read_maps}"
            fi

            cat !{isoform_fastas} > collapse_combined.fasta
            cat !{isoform_beds} > collapse_combined.bed
            cat !{isoform_gtfs} > collapse_combined.gtf
            cat !{isoform_read_maps} > collapse_combined.read.map.txt
            
        '''
}