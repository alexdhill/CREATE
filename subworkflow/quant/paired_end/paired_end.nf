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
import java.nio.file.Files
import java.nio.file.Paths


include { prefetch } from "../../../modules/sra/prefetch/prefetch.nf"
include { fasterq_dump_paired } from "../../../modules/sra/fasterq-dump/fasterq-dump_paired.nf"
include { count_reads_pe } from "../../../modules/bash/count_reads/count_reads_pe.nf"
include { trim_reads_paired } from "../../../modules/trim-galore/trim_reads/trim_reads_paired.nf"
include { salmon_quant_paired } from "../../../modules/salmon/salmon_quant/salmon_quant_paired.nf"
include { compile_quantifications } from "../../../modules/R/compile_quantifications/compile_quantifications.nf"
include { run_analysis } from "../../../modules/R/run_analysis/run_analysis.nf"

workflow PAIRED_END
{
    take:
        is_acc
    main:
        if (params.library=="paired_end")
        {
            reference = Channel.fromPath(params.ref)
            metadata = Channel.fromPath(params.metadata)
            if (is_acc)
            {
                log.info("Downloading reads before running...")
                Channel.from(file(params.samples).readLines())
                | prefetch
                | fasterq_dump_paired
                | set{reads}
            } else
            {
                log.info("Parsing reads from file...")
                reads = Channel.fromFilePairs(params.samples+"/"+params.pattern)
                    .map{sample -> [sample[0], sample[1][0], sample[1][1]]}
            }

            reads
            | count_reads_pe
            | trim_reads_paired
            | combine(reference)
            | salmon_quant_paired
            | collect
            | map{quants -> [quants]}
            | combine(reference)
            | compile_quantifications
            | combine(metadata)
            | run_analysis
        }
}