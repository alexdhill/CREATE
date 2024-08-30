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


/*
 * Add pattern import to compile the file matcher
 */
import java.nio.file.FileSystems
import java.io.File
import java.util.Arrays

include { PAIRED_END } from '../../subworkflow/quant/paired_end/paired_end.nf'
include { SINGLE_END } from '../../subworkflow/quant/single_end/single_end.nf'
include { SINGLE_CELL } from '../../subworkflow/quant/single_cell/single_cell.nf'
include { NANOPORE } from '../../subworkflow/quant/nanopore/nanopore.nf'

workflow QUANT
{
    /*
     * File matcher for debug
     */
    sampleDir = (params.samples.lastIndexOf("/")+1==params.samples.length())?params.samples:params.samples+"/"
    log.info("Listing ${sampleDir}${params.pattern}")
    glob = FileSystems.getDefault().getPathMatcher("glob:${sampleDir}${params.pattern}")
    files = Arrays.stream(new File(sampleDir).listFiles())
        .filter(f -> glob.matches(f.toPath()))
        .toArray();
    Arrays.sort(files)
    pair = (["single_end", "nanopore"].contains(params.library))?1:2
    Arrays.copyOfRange(files, 0, Math.min(files.length, 10)).eachWithIndex{f,i -> {
        if (i%(pair*2)<pair) log.info(">   \u001B[36m${f.toString()}\u001B[0m")
        else log.info(">   ${f.toString()}")
    }}
    if (files.size() > 10) log.info(">   . . . and \u001B[32m${files.size()-10}\u001B[0m more")
    log.info(" ")


    /*
     * Gather files as channels
     */
    if (["paired_end","single_cell"].contains(params.library))
    {
        reads = Channel.fromFilePairs(params.samples+"/"+params.pattern)
            .map{sample -> [sample[0], sample[1][0], sample[1][1]]}
    }
    else
    {
        reads = Channel.fromPath(params.samples+"/"+params.pattern)
            .map{sample -> [sample.name.split(/\.f(ast)?q(\.gz)?/)[0], sample]}
    }

    /*
     * Run subworkflows
     */
    switch (params.library)
    {
        case "paired_end":
            PAIRED_END(reads)
        case "single_end":
            SINGLE_END(reads)
        case "single_cell":
            SINGLE_CELL(reads)
        case "nanopore":
            NANOPORE(reads)
    }
}