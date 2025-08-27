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


include { prefetch } from "../../../modules/sra/prefetch/prefetch.nf"
include { fasterq_dump_single } from "../../../modules/sra/fasterq-dump/fasterq-dump_single.nf"
include { count_reads_se } from "../../../modules/bash/count_reads/count_reads_se.nf"
include { trim_reads_single } from "../../../modules/trim-galore/trim_reads/trim_reads_single.nf"
include { salmon_quant_single } from "../../../modules/salmon/salmon_quant/salmon_quant_single.nf"
include { compile_quantifications } from "../../../modules/R/compile_quantifications/compile_quantifications.nf"
include { run_analysis } from "../../../modules/R/run_analysis/run_analysis.nf"

workflow SINGLE_END
{
    take:
        is_acc
    main:
        if (params.library=="single_end")
        {
            reference = Channel.fromPath(params.ref)
            parameters = Channel.fromPath(params.parameters)
            metadata = Channel.fromPath(params.metadata)
            if (is_acc)
            {
                log.info("Downloading reads before running...")
                Channel.fromPath(file(params.samples).readLines())
                | prefetch
                | fasterq_dump_single
                | set{reads}
            } else
            {
                reads = Channel.fromPath(params.samples+"/"+params.pattern)
                .map{sample -> [sample.name.split(/\.f(ast)?q(\.gz)?/)[0], sample]}
            }

            reads
            | count_reads_se
            | concat(reads)
            | groupTuple()
            | map{ res -> [res[0], res[1][0], res[1][1]] }
            | combine(parameters)
            | trim_reads_single
            | combine(reference)
            | combine(parameters)
            | salmon_quant_single
            | collect
            | map{quants -> [quants]}
            | combine(reference)
            | combine(metadata)
            | compile_quantifications
            | run_analysis
        }
}