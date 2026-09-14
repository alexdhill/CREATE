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


include { nanoplot_report } from "../../../modules/nanoplot/nanoplot_report/nanoplot_report.nf"
include { multiqc_report_long } from "../../../modules/multiqc/multiqc_report/multiqc_report_long.nf"

workflow NANOPORE
{
    take:
        raw_reads
        trimmed_reads
        align_logs
        quants
    main:
        if (params.parameters) parameters = Channel.fromPath(params.parameters)
        else parameters = Channel.fromPath(projectDir + "/assets/NULL")

        raw_reads
        | map{sample -> [sample[0], "raw", sample[1]]}
        | concat(trimmed_reads.map{sample -> [sample[0], "trimmed", sample[1]]})
        | combine(parameters)
        | nanoplot_report
        | collect
        | map{reports -> [reports]}
        | combine(align_logs.collect().map{logs -> [logs]})
        | combine(quants)
        | combine(parameters)
        | multiqc_report_long
}
