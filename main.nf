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
 * Imports
 */
import groovy.util.logging.Slf4j
// ERROR -> WARN -> DEBUG -> INFO


/*
 * Handle parameters
 */

/*
 * Create parameter set
 */
params
{
    help = false
}
if (params.help)
{
    println '''Usage: nextflow run alexdhill/CREATE [--quant/--reference/--analyze] [options]
    reference         Build a reference for CREATE
    quant             Quantify RNA-seq reads with CREATE

options:
    Globals:
        --help        Print this help message
        --container   Containerization method [docker(default),conda] TODO!
        --outdir      Output directory [REQUIRED]
        --log         Log level [OFF(default),INFO,DEBUG]
        --exec        Executor mode [local(default),slurm]
        --force       Overwrite existing files [default=false]
        --keep        Keep intermediate files [default=false]
        --limits      Add sample-based resources limits when running locally [default=false]
        --clean       Clean up working directory (breaks `-resume`!) [default=false]

    --reference
        --genome      Genome assembly [HG38,T2T] [REQUIRED]
        --version     Gencode version [HG38 only] [default=39]
        --index       Desired index [short,long,single_cell] [default=short]
        --isoquant    Novel isoform annotation [BETA]

    --quant [default]
        --samples     Input samples directory [REQUIRED]
        --ref         CREATE reference directory [REQUIRED]
        --library     Library type [paired_end,nanopore,single_cell] [REQUIRED]
                      NOTE: the reference provided must contain a compatible index
        --pattern     File pattern for input samples
        --metadata    An unnamed, two-column CSV file of samples,condition

    Subworkflow-specific options:
        [--exec local]
            --threads     Maximum number of threads to allocate TODO!
            --memory      Maximum memory to allocate TODO!
        [--exec slurm]
            --njobs       Maximum number of jobs to run concurrently [default=15]
            --scratch     Scratch directory for temporary files [default=false]
            --account     The user account to submit jobs
        [--quant --library single_cell]
            --barcodes    A list (txt) of all barcodes used in the experiment [REQUIRED]
            --chemistry   The single-cell chemsitry used [dropseq,chromium,chromiumV3] [default=chromiumV3]
        [--quant --library nanopore]
            --dcs         File with the DCS/RCS sequence used (for removal) [REQUIRED]
            --discovery   Discover novel isoforms [default=false] [BETA]
'''
    System.exit 0
}
params
{
    container = 'docker'
    outdir = ''
    log = 'OFF'
    exec = 'local'
    clean = false
    keep = false
    limits = false
    force = false
    quant = true
    reference = false
}
if (params.reference)
{
    params
    {
        genome = ''
        isoquant = false
        index = ''
        quant = false
        params.version = ''
    }
    if (params.genome=='HG38')
    {
        params
        {
            version = '39'
        }
    }
    else if (params.genome!="T2T")
    {
        println("--genome: ${params.genome}");
        throw new IllegalArgumentException('Genome must be either "HG38" or "T2T"')
    }
    if (params.genome=="T2T" && params.version!='')
    {
        println("--version: ${params.version}");
        throw new IllegalArgumentException('Gencode version cannot be specified for T2T genome')
    }
    index = new ArrayList(Arrays.asList(params.index.split(',')))
    index.retainAll(['short','long','single_cell'])
    if (index.size()!=params.index.split(',').size())
    {
        println("--index: ${params.index}");
        throw new IllegalArgumentException('index must be one or more of: short, long, single_cell')
    }
}
else if (params.quant)
{
    params
    {
        samples = ''
        ref = ''
        library = 'paired_end'
        metadata = ''
    }
    if (['paired_end','single_cell'].contains(params.library))
    {
        params.pattern = '*_R{1,2}_*.fastq.gz'
        if (params.library=="single_cell")
        {
            params
            {
                barcodes = ''
                chemistry = "chromiumV3"
            }
            if (params.barcodes=='')
            {
                println("--barcodes: ${params.barcodes}");
                throw new IllegalArgumentException('Barcode list must be specified')
            }
            if (!["chromiumV3","chromium","dropseq"].contains(params.chemistry))
            {
                println("--chemistry: ${params.chemistry}");
                throw new IllegalArgumentException('Chemistry must be one of: chromiumV3, chromium, dropseq')
            }
        }
    }
    else if (params.library=='nanopore')
    {
        params.pattern = '*.fastq.gz'
        params
        {
            dcs = ""
            discovery = false
        }
        if (params.dcs=='')
        {
            println("--dcs: ${params.dcs}");
            throw new IllegalArgumentException('DCS/RCS sequence file must be specified')
        }
    }
    else if (params.library=='single_end')
    {
        params.pattern = '*_R1_*.fastq.gz'
    }
    else
    {
        println("--library: ${params.library}");
        throw new IllegalArgumentException('Library type must be either "paired_end", "single_end", "nanopore", or "single_cell"')
    }
    if (params.samples=='')
    {
        println("--samples: ${params.samples}");
        throw new IllegalArgumentException('Input sample directory must be specified')
    }
    if (params.metadata=='')
    {
        println("--metadata: ${params.metadata}");
        throw new IllegalArgumentException('Sample metadata file must be specified')
    }
    if (params.ref=='')
    {
        println("--ref: ${params.ref}");
        throw new IllegalArgumentException('CREATE reference directory must be specified')
    }
}
else
{
    throw new IllegalArgumentException('CREATE-seq mode must be specified either "--quant" or "--reference"')
}
if (params.reference && params.quant)
{
    println("--quant: ${params.quant}");
    println("--reference: ${params.reference}");
    throw new IllegalArgumentException('CREATE-seq mode must be specified either "--quant" XOR "--reference"')
}
if (!['OFF','INFO','DEBUG'].contains(params.log))
{
    println("--log: ${params.log}");
    throw new IllegalArgumentException('Log level must be one of: OFF, INFO, DEBUG')
}
if (params.outdir=='')
{
    println("--outdir: ${params.outdir}");
    throw new IllegalArgumentException('Output directory must be specified')
}
if (params.exec=="local")
{
    if (params..get('scratch'))
    {
        throw new IllegalArgumentException('Cannot assign scratch directory in local mode')
    }
    if (params.get('account'))
    {
        throw new IllegalArgumentException('Cannot assign account in local mode')
    }
    if (params.get('njobs'))
    {
        throw new IllegalArgumentException('Cannot assign number of jobs in local mode')
    }
}
else if (params.exec=="slurm")
{
    if (params.get('threads'))
    {
        throw new IllegalArgumentException('Cannot assign threads in slurm mode')
    }
    if (params.get('memory'))
    {
        throw new IllegalArgumentException('Cannot assign memory in slurm mode')
    }
}


/*
 * Apply parameters
 */
process {
    cache = true
    errorStrategy = 'finish'
    beforeScript = 'chmod 755 .'
}

cleanup = params.clean
if (params.container=='conda')
{
    conda.enabled = true
    process.conda = '${projectDir}/bin/conda/create_env.yaml'
    params.docker = false
}
else if (params.container=='docker')
{
    docker.enabled = true
    docker.runOptions = '-u \$(id -u):\$(id -g)'
    process.container = 'alexdhill/create:test'
}
else
{
    throw new IllegalArgumentException('Container must be either "docker" or "conda"')
}
if (params.exec=='slurm')
{
    params {
        njobs = 15
        scratch = false
        account = ''
    }
    process.executor = 'slurm'
    process.scratch = params.scratch
    executor.queueSize = params.njobs
    if (params.account!='') executor.account = params.account
}
else if (params.exec=="local")
{
    if (params.get('threads')==null || params.threads<0)
    {
        println("--threads: ${params.threads}");
        throw new IllegalArgumentException('Maximum number of threads must be at least 0')
    }
    else if (params.threads>0)
    {
        executor.cpus = params.threads
    }
    if (params.get('memory')==null || params.memory<0)
    {
        println("--memory: ${params.memory}");
        throw new IllegalArgumentException('Memory must be at least 0 GB')
    }
    else if (params.memory>0)
    {
        executor.memory = params.memory
    }
}
else
{
    println("--exec: ${params.exec}");
    throw new IllegalArgumentException('Execution mode must be either "local" or "slurm"')
}

/*
 * Print final parameters
 */
def print_val(k,v)
{
    type=v.getClass().getSimpleName()
    key = k
    while (key.length()<12) key += " "
    switch(type)
    {
        case "Boolean":
            log.info("    $key => \u001B[32m$v\u001B[0m")
            return
        case "String":
            value = v
            value = value.replace(',', '\u001B[0m,\u001B[33m')
            value = value.replace('*', '\u001B[35m*\u001B[33m')
            log.info("    $key => \u001B[33m$value\u001B[0m")
            return
        case "Integer":
            log.info("    $key => \u001B[34m$v\u001B[0m")
            return
        case "Float":
            log.info("    $key => \u001B[34m$v\u001B[0m")
            return
        default:
            log.info("    $key => $v")
            return
    }
}
log.info("Run parameters:")
params.each{k, v -> ['',false,0].contains(v)?:print_val(k,v)}
log.info(" ")
params.manage_resources = (params.limits || params.exec!='local')

/*
 * Handle post-print params
 */
if (params.get('reference'))
{
    if (![null,false].contains(params.get('isoquant'))) params.version = '00'
    else if (params.get('genome')=='T2T') params.version = '2'
}


/*
 * Make CREATE workflow selector
 */
include { REFERENCE } from './workflow/reference/reference.nf'
include { QUANT } from './workflow/quant/quant.nf'

workflow CREATE
{
    if (params.get('quant')!=null && params.get('quant')==true) QUANT()
    else REFERENCE()
}


/*
 * Run CREATE
 */
workflow { CREATE() }