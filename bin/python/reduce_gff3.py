#!/bin/python3

from sys import stdin

chrOrder = {
    "chr1": 0,
    "chr2": 1,
    "chr3": 2,
    "chr4": 3,
    "chr5": 4,
    "chr6": 5,
    "chr7": 6,
    "chr8": 7,
    "chr9": 8,
    "chr10": 9,
    "chr11": 10,
    "chr12": 11,
    "chr13": 12,
    "chr14": 13,
    "chr15": 14,
    "chr16": 15,
    "chr17": 16,
    "chr18": 17,
    "chr19": 18,
    "chr20": 19,
    "chr21": 20,
    "chr22": 21,
    "chrX": 22,
    "chrY": 23,
    "chrM": 24,
}


class Exon:
    def __init__(self, start, stop):
        self.start = start
        self.stop = stop


class Transcript:
    def __init__(self, id, src, type, start, stop):
        self.id = id
        self.src = src
        self.type = type
        self.start = start
        self.stop = stop
        self.exons = []


class Gene:
    def __init__(self, id, src, name, type, chr, build, start, stop, strand):
        self.id = id
        self.src = src
        self.name = name
        self.type = type
        self.chr = chr
        self.build = build
        self.start = start
        self.stop = stop
        self.strand = strand
        self.transcripts = {}


def reduce_gff3():
    genes = {}
    for line in stdin:
        if line.startswith("#"):
            continue
        (chr, build, feat, start, stop, _, strand, _, attr) = line.strip().split("\t")
        if feat == "gene":
            attrs = attr.split(";")
            gid, gname, gtype, srcid = "", "", "", ""
            for attr in attrs:
                if attr.strip().startswith("gene_id"):
                    gid = attr.split()[1].replace('"', "")
                elif attr.strip().startswith("gene_name"):
                    gname = attr.split()[1].replace('"', "")
                elif attr.strip().startswith("gene_biotype"):
                    gtype = attr.split()[1].replace('"', "")
                elif attr.strip().startswith("source_gene"):
                    srcid = attr.split()[1].replace('"', "")
            genes[gid] = Gene(
                id=gid,
                src=srcid,
                name=gname,
                type=gtype,
                chr=chr,
                build=build,
                start=int(start),
                stop=int(stop),
                strand=strand,
            )
        elif feat == "transcript":
            attrs = attr.split(";")
            gid, txid, txtype, srctx = "", "", "", ""
            for attr in attrs:
                if attr.strip().startswith("gene_id"):
                    gid = attr.split()[1].replace('"', "")
                elif attr.strip().startswith("transcript_id"):
                    txid = attr.split()[1].replace('"', "")
                elif attr.strip().startswith("transcript_biotype"):
                    txtype = attr.split()[1].replace('"', "")
                elif attr.strip().startswith("source_transcript"):
                    srctx = attr.split()[1].replace('"', "")
            genes[gid].transcripts[txid] = Transcript(
                id=txid, type=txtype, src=srctx, start=int(start), stop=int(stop)
            )
        elif feat == "exon":
            attrs = attr.split(";")
            gid, txid = "", ""
            for attr in attrs:
                if attr.strip().startswith("gene_id"):
                    gid = attr.split()[1].replace('"', "")
                elif attr.strip().startswith("transcript_id"):
                    txid = attr.split()[1].replace('"', "")
            genes[gid].transcripts[txid].exons.append(
                Exon(start=int(start), stop=int(stop))
            )
    for gene in sorted(genes.values(), key=lambda g: (chrOrder[g.chr], g.start)):
        print(
            "\t".join(
                [
                    gene.chr,
                    gene.build,
                    "gene",
                    str(gene.start),
                    str(gene.stop),
                    ".",
                    gene.strand,
                    ".",
                    'gene_id "{}"; gene_name "{}"; gene_biotype "{}"; source_gene "{}"; '.format(
                        gene.id, gene.name, gene.type, gene.src
                    ),
                ]
            )
        )
        for tx in sorted(gene.transcripts.values(), key=lambda t: t.start):
            print(
                "\t".join(
                    [
                        gene.chr,
                        gene.build,
                        "transcript",
                        str(tx.start),
                        str(tx.stop),
                        ".",
                        gene.strand,
                        ".",
                        'gene_id "{}"; gene_name "{}"; gene_biotype "{}"; source_gene "{}"; transcript_id "{}"; transcript_biotype "{}"; source_transcript "{}"; '.format(
                            gene.id,
                            gene.name,
                            gene.type,
                            gene.src,
                            tx.id,
                            tx.type,
                            tx.src,
                        ),
                    ]
                )
            )
            for exon in sorted(tx.exons, key=lambda e: e.start):
                print(
                    "\t".join(
                        [
                            gene.chr,
                            gene.build,
                            "exon",
                            str(exon.start),
                            str(exon.stop),
                            ".",
                            gene.strand,
                            ".",
                            'gene_id "{}"; gene_name "{}"; gene_biotype "{}"; source_gene "{}"; transcript_id "{}"; transcript_biotype "{}"; source_transcript "{}"; '.format(
                                gene.id,
                                gene.name,
                                gene.type,
                                gene.src,
                                tx.id,
                                tx.type,
                                tx.src,
                            ),
                        ]
                    )
                )


if __name__ == "__main__":
    reduce_gff3()
