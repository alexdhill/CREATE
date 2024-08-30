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


include { make_splintr_transcripts } from '../../../modules/R/make_splintr_transcripts/make_splintr_transcripts.nf'
include { salmon_index } from '../../../modules/salmon/salmon_index/salmon_index.nf'
include { link_transcriptome } from '../../../modules/R/link_transcriptome/link_transcriptome.nf'
include { make_database } from '../../../modules/R/make_database/make_database.nf'

workflow SINGLE_CELL
{
    take:
        complete // [complete.gtf, complete.fa, genome]
    main:
        if (params.index.split(',').contains('single_cell'))
        {
            complete
            | map{dat -> [dat[0], dat[2]]}
            | make_splintr_transcripts
            | map{dat -> dat[0]}
            | combine(["splintr"])
            | salmon_index
            
            // complete
            // | map{dat -> [dat[1], dat[2]]}
            // | combine(salmon_index.out)
            // | link_transcriptome & make_database
        }
}