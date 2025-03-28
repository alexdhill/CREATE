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
import groovy.util.logging.Slf4j

include { download_reference } from '../../modules/bash/download_reference/download_reference.nf'
include { download_gencode_transcripts } from '../../modules/bash/download_gencode_transcripts/download_gencode_transcripts.nf'
include { download_gencode_annotation } from '../../modules/bash/download_gencode_annotation/download_gencode_annotation.nf'
include { download_repeat_regions } from '../../modules/bash/download_repeat_regions/download_repeat_regions.nf'
include { make_repeat_transcripts } from '../../modules/bedtools/make_repeat_transcripts/make_repeat_transcripts.nf'
include { download_repeat_annotation } from '../../modules/R/download_repeat_annotation/download_repeat_annotation.nf'
include { make_complete_transcripts } from '../../modules/bash/make_complete_transcripts/make_complete_transcripts.nf'
include { make_complete_annotation } from '../../modules/bash/make_complete_annotation/make_complete_annotation.nf'
include { make_transcript_map } from '../../modules/R/make_transcript_map/make_transcript_map.nf'
include { link_transcriptome } from '../../modules/R/link_transcriptome/link_transcriptome.nf'
include { make_database } from '../../modules/R/make_database/make_database.nf'

include { SHORT } from '../../subworkflow/reference/short/short.nf'
include { LONG } from '../../subworkflow/reference/long/long.nf'
include { SINGLE_CELL } from '../../subworkflow/reference/single_cell/single_cell.nf'
include { DISCOVER } from '../../subworkflow/reference/discover/discover.nf'

workflow REFERENCE
{
    download_reference() // get reference genome
    download_gencode_annotation()
    download_repeat_regions()
    if (params.genome!="T2T")
    {
        log.info("Non-T2T genome detected. Downloading txome")
        Channel.from(["null.fa", "null.crt"])
        | download_gencode_transcripts
    } else {
        log.info("T2T genome detected. Generating txome")
        download_gencode_annotation.out
        | combine(download_reference.out)
        | download_gencode_transcripts
    }

    log.info("Starting main logic...")
    download_reference.out
    | combine(download_repeat_regions.out) // [genome.fa, repeats.bed]
    | (download_repeat_annotation & make_repeat_transcripts)
    | concat // repeats.gtf, repeats.fa
    | collect // [repeats.gtf, repeats.fa]
    | combine(download_gencode_annotation.out) // [repeats.gtf, repeats.fa, gencode.gtf]
    | combine(download_gencode_transcripts.out) // [repeats.gtf, repeats.fa, gencode.gtf, gencode.fa]
    | (make_complete_annotation & make_complete_transcripts)
    | concat // complete.gtf, complete.fa
    | collect // [complete.gtf, complete.fa]
    | combine(download_reference.out) // [complete.gtf, complete.fa, genome.fa]
    | link_transcriptome & SHORT & LONG & SINGLE_CELL

    download_reference.out
    | DISCOVER

    download_gencode_annotation.out
    | combine(download_repeat_annotation.out)
    | make_transcript_map
}
