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
include { generate_transcriptome } from '../../modules/bash/generate_transcriptome/generate_transcriptome.nf'
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

NOENT = projectDir+'/assets/NULL'

workflow REFERENCE
{
    // Get reference genome
    download_reference()
    | set{reference}

    // Get RepeatMasker track
    download_repeat_regions()
    | set{repeatmasker}

    // Get gene annotation and transcriptome
    if (params.isoquant) {
        log.info("Using provided annotation")
        reference
        | combine(Channel.fromPath(params.isoquant))
        | generate_transcriptome
        | set{transcriptome}

        annotation = Channel.fromPath(params.isoquant)
    } else {
        log.info("Getting annotation")
        reference
        | combine(download_gencode_annotation())
        | download_gencode_transcripts
        | set{transcriptome}

        download_gencode_annotation.out
        | set{annotation}
    }

    log.info("Starting main logic...")
    reference
    | combine(repeatmasker) // [genome.fa, repeats.bed]
    | (download_repeat_annotation & make_repeat_transcripts)
    | concat // repeats.gtf, repeats.fa
    | collect // [repeats.gtf, repeats.fa]
    | combine(annotation) // [repeats.gtf, repeats.fa, gencode.gtf]
    | combine(transcriptome) // [repeats.gtf, repeats.fa, gencode.gtf, gencode.fa]
    | (make_complete_annotation & make_complete_transcripts)
    | concat // complete.gtf, complete.fa
    | collect // [complete.gtf, complete.fa]
    | combine(reference) // [complete.gtf, complete.fa, genome.fa]
    | link_transcriptome & SHORT & LONG & SINGLE_CELL

    reference
    | DISCOVER

    annotation
    | combine(download_repeat_annotation.out)
    | make_transcript_map
}