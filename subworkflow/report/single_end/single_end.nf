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


include { fastqc_report_single } from "../../../modules/fastqc/fastqc_report/fastqc_report_single.nf"
include { multiqc_report_short } from "../../../modules/multiqc/multiqc_report/multiqc_report_short.nf"

workflow SINGLE_END
{
    take:
        raw_reads
        trimmed_reads
        quants
    main:
        if (params.parameters) parameters = Channel.fromPath(params.parameters)
        else parameters = Channel.fromPath(projectDir + "/assets/NULL")

        raw_reads
        | map{sample -> [sample[0], "raw", sample[1]]}
        | concat(trimmed_reads.map{sample -> [sample[0], "trimmed", sample[1]]})
        | combine(parameters)
        | fastqc_report_single
        | collect
        | map{reports -> [reports]}
        | combine(quants)
        | combine(parameters)
        | multiqc_report_short
}
