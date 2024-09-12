process
{
    // Ref

    withName: "download_gencode_annotation"
    {
        queue = "short"
        time = '5m'
    }
    withName: "download_gencode_transcripts"
    {
        queue = "short"
        time = '5m'
    }
    withName: "download_reference"
    {
        queue = "short"
        time = '10m'
    }
    withName: "download_repeat_regions"
    {
        queue = "short"
        time = '5m'
    }

    withName: "download_repeat_annotation"
    {
        queue = "short"
        time = '45m'
    }

    withName: "make_complete_annotation"
    {
        queue = "short"
        time = '5m'
    }
    withName: "make_complete_transcripts"
    {
        queue = "short"
        time = '5m'
    }

    withName: "link_transcriptome"
    {
        queue = "medium"
        time = '90m'
    }

    withName: "star_index_genome"
    {
        queue = "long"
        time = '24h'
    }

    withName: "make_splintr_transcripts"
    {
        queue = "medium"
        time = '3h'
    }

    withName: "salmon_index"
    {
        queue = "short"
        time = '1h'
    }

    withName: "minimap2_index"
    {
        queue = "short"
        time = '5m'
    }

    // Quant

    withName: "count_reads_pe"
    {
        queue = "short"
        time = '10m'
    }
    withName: "count_reads_se"
    {
        queue = "short"
        time = '10m'
    }
    withName: "count_reads_np"
    {
        queue = "short"
        time = '10m'
    }

    withName: "trim_reads_paired"
    {
        queue = "short"
        time = '1h'
    }
    withName: "trim_reads_single"
    {
        queue = "short"
        time = '1h'
    }
    withName: "trim_reads_nanopore"
    {
        queue = "short"
        time = '1h'
    }

    withName: "salmon_quant_paired"
    {
        queue = "medium"
        time = '6h'
    }
    withName: "salmon_quant_single"
    {
        queue = "medium"
        time = '6h'
    }
    withName: "salmon_quant_nanopore"
    {
        queue = "medium"
        time = '6h'
    }

    withName: "compile_quantifications"
        {
            queue = "medium"
            time = '3h'
        }
}
