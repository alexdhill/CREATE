<!--
  REQUIRED NOTICE: Copyright (c) 2020-2023, Regents of the University of California
  All rights reserved. https://polyformproject.org/licenses/noncommercial/1.0.0
  
  This software was developed by the Daniel Kim lab at the University of California, Santa Cruz.
  Authors: Roman E. Reggiardo, Vikas Peddu, Alex D. Hill
  
  The licensor grants you a copyright license for the software to do everything you might do with
  the software that would otherwise infringe the licensor’s copyright in it for any permitted
  purpose.
  
  As far as the law allows, the software comes as is, without any warranty or condition, and the
  licensor will not be liable to you for any damages arising out of these terms or the use or
  nature of the software, under any kind of legal claim.
-->


# `CREATE`-seq
## Comprehensive REptitive And Transposable Element Sequencing
### Alex D. Hill

CREATE is a pipeline for the quantification of repetitive elements alongside canonical genes. It is designed to be very high throughput, and can be run very quickly on a large number of samples. The core of CREATE's functions are based on Trim-Galore, Salmon, and Tximeta for fast read processing and storage.

There are three modes to CREATE:
  - REFERENCE: This mode creates a rich 'complete' set of reference files that contain a formatted compilation of both GENCODE annotations and transcripts, as well as repetitive element annotations and transcripts pulled from the RepeatMasker track on the UCSC genome browser. Reference directories can also be downloaded for ease of use.
  - QUANT: This mode takes a set of raw reads and outputs a highly efficient H5+RDS file that efficiently store a SummarizedExperimet (or SingleCellExperiment) object that can be used for downstream analysis.
  - DISCOVER: This mode will take an existing CREATE reference and use long read samples to create a enriched reference with novel isoforms, and then quantify short read samples with the novel reference.

#### Installation:
As a nextflow pipeline, CREATE can be run with only an installation of conda. Make an environment with Nextflow installed and force `CREATE` to utilize conda and the container format.  
` $ conda create -n CREATE -c bioconda nextflow`  
The pipeline can then be run by calling: ` $ nextflow run alexdhill/CREATE`  
*NOTE*: a docker installation is highly recommended (any used by default) to guarantee compatibility with most systems. 

#### usage:
<!--    discover          Discover novel isoforms [BETA] -->
```
 $ nextflow run alexdhill/CREATE --help
Usage: nextflow run alexdhill/CREATE [--quant/--reference/--discover] [options]
    reference         Build a reference for CREATE
    quant             Quantify RNA-seq reads with CREATE
    discover          Discover novel isoforms [BETA]

options:
    Globals (used in all modes):
        --help        Print this help message
        --container   Containerization method [docker(default),conda]
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

    --quant
        --samples     Input samples directory [REQUIRED]
        --ref         CREATE reference directory [REQUIRED]
        --library     Library type [paired_end,nanopore,single_cell] [REQUIRED]
                      NOTE: the reference provided must contain a compatible index
        --pattern     File pattern for input samples
        --metadata    An unnamed, two-column CSV file of samples,condition

    --discover
        --ref          CREATE reference directory [REQUIRED]
        --dcs          The DNA/RNA control sequence fasta file [REQUIRED]
        --prefixes     A list of file prefixes for the paired long-short reads [REQUIRED]
        --long_reads   The directory of the long-read sequences [REQUIRED]
        --paired_reads The directory of the paired-end sequences [REQUIRED]
        --lr_pattern   The pattern for the long-read files [default='*.fastq.gz']
        --pe_pattern   The pattern for the paired-end files [default='*_R{1,2}_*.fastq.gz']
        --dump         Save a 'novel' CREATE reference [default=false]

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
            --chemistry   The single-cell chemistry used [dropseq,chromium,chromiumV3(default)]
        [--quant --library nanopore]
            --dcs         File with the DCS/RCS sequence used (for removal) [REQUIRED]
```

Global arguments are used in all modes to control how nextflow runs `CREATE`. For each individual mode, there are certain arguments that are reuquired for input data. Nested arguments shown above show what arguments are allowed (or required) depending on the running mode.

#### examples: REFERENCE
To build a new `CREATE` reference directory for illumina quantification, you can run:
```
nextflow run alexdhill/CREATE --reference \
  --genome HG38 \
  --version 39 \
  --index short \
  --outdir HG38_reference.crt
```

#### examples: QUANT
```
nextflow run alexdhill/CREATE --quant \
  --samples <path/to/samples> \
  --ref <path/to/reference> \
  --library paired_end \
  --outdir create_results
```

#### examples: DISCOVER
```
nextflow run alexdhill/CREATE --discover \
  --long_reads <path/to/nanopore> \
  --paired_reads <path/to/illumina> \
  --prefixes <prefixes.txt> \
  --dcs <dcs_sequence.fa> \
  --ref <path/to/reference>
```

#### Flags and arguments
 - `--account`: specify the account username. [ONLY used with `--exec slurm`]
 - `--barcodes`: path to a newline-delimited file of barcode sequences [ONLY used with `--quant` **and** `--library single_cell`]
 - `--chemistry`: the single-cell chemstry to use [select from chromium, chromiumV3, or dropseq] [ONLY used with `--library single_cell`]
 - `--clean`: removes the work directory upon completion of workflow [BETA]
 - `--container`: type of environment manager to use [select from docker or conda]
 - `--dcs`: path to a fasta file with the DCS or RCS sequence used [ONLY used with `--quant` **and** `--library nanopore`]
 - `--discover`: specifies to run the DISCOVER pipeline to generate a novel isoform enriched CREATE directory and quantify reads using matched nanopore and paired-end sequencing data
 - `--dump`: path to save the novel CREATE reference generated from the discover pipeline [ONLY used with `--discover`]
 - `--exec`: specifies the method of job execution [select from local or slurm]
 - `--force`: overwrite existing output if possible
 - `--genome`: specifies which reference genome to use [select from HG38, MM39, or T2T] [ONLY used with `--reference`]
 - `--help`: prints help message and exits
 - `--index`: specifies which quantification index(es) to generate. Multiple can be specified as comma-delimited [select from short, long, or single_cell] [ONLY used with `--reference`]
 - `--isoquant`: path to an isoquant-generated GTF to use instead of the specified reference genome's GENCODE GTF. [ONLY used with `--reference`]
 - `--keep`: keeps generated intermediate files
 - `--library`: specifies the type of cDNA library generated [select from paired_end, single_end, single_cell, or nanopore] [ONLY used with `--quant`]
 - `--limits`: specifies whether to limit the resources provided for each job [BETA] [ONLY used with `--exec local`]
 - `--log`: specifies the logging level. This will print the arguments and found samples to the console, and save each individual job's logging to the ".command.log" file in the work directory.
 - `--long_reads`: specifies the directory to search for the nanopore read samples [ONLY used with `--discover`]
 - `--lr_pattern`: specifies a glob pattern to select the desired nanopore read files [ONLY used with `--discover`]
 - `--memory`: specifies the mamimum memory that can be allocated by **all** active jobs together [BETA]
 - `--metadata`: path to a comma-separated file matching sample names and their condition in that order. Sample names are determined by the pattern used to select the files. For paired reads, the default uses the full sample name *until* _R{1,2} [ONLY used with `--quant`]
 - `--njobs`: the maximum number of jobs to spawn in parallel [ONLY used with `--exec slurm`]
 - `--outdir`: path to the output directory to save all results
 - `--pattern`: the pattern used to match the desired samples from the `--samples` directory [ONLY used with `--quant`]
 - `--paired_reads`: the directory to search for the paired-end samples [ONLY used with `--discover`]
 - `--pe_pattern`: the glob pattern to use to select the desired reads from the paired-end reads directory [ONLY used with `--discover`]
 - `--prefixes`: path to a newline-delimited file of all the matched nanopore and paired-end sample prefixes [ONLY used with `--discover`]
 - `--quant`: run the QUANT pipeline to quantify RNA-seq samples with a CREATE reference
 - `--reference`: run the REFERENCE pipeline to build a CREATE reference directory
 - `--ref`: path to the CREATE reference directory to use [ONLY used with `--quant` or `--discover`]
 - `--samples`: path to the directory to search for the RNA-seq samples [ONLY used with `--quant`]
 - `--scratch`: path to the scratch directory to store working data in [ONLY used with `--exec slurm`]
 - `--threads`: the maximum number of threads to allocate between **all** running jobs [BETA]
 - `--version`: the version of the gencode annotation to use, must be prefixed by `M` if building a mouse reference [ONLY used with `--reference` and `--genome HG38`]

#### Roadmap:
- [X] Finish Reference workflows
  - [X] Finish salmon reference for short reads
  - [X] Finish minimap reference for long reads
  - [X] Finish splintr reference for single-cell
- [X] Finish Quant workflows
  - [X] Finish Paired read quant
  - [X] Finish Single read quant
  - [X] Finish Single-Cell quant
  - [X] Finish Nanopore quant
    - [X] Add DCS/RCS removal
- [X] Add novel isoform discovery (FLAIR)
- [X] Add Analysis ~~workflow~~ step to quant
  - [X] Finish biotype plots
  - [X] Finish PCA plots
  - [X] Finish Volcano plots
- [X] Flesh out README
- [ ] ~~Add profile selection to nextflow config~~ Add read-based resource requirements to processes
- [ ] Add accession list support to automate downloads [feature-samplesheet]
- [ ] Add semi-intelligent prefix parsing for compilation [feature-metadata]
- [ ] Add variance-stabilized normalized counts to output assays [feature-metadata]
- [ ] Add alignment mode [feature-align]
- [ ] Add mouse genome (GRCm39/mm39) support [feature-mouse]