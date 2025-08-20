#!/bin/python3

chrOrder = {
    'chr1': 0,
    'chr2': 1,
    'chr3': 2,
    'chr4': 3,
    'chr5': 4,
    'chr6': 5,
    'chr7': 6,
    'chr8': 7,
    'chr9': 8,
    'chr10': 9,
    'chr11': 10,
    'chr12': 11,
    'chr13': 12,
    'chr14': 13,
    'chr15': 14,
    'chr16': 15,
    'chr17': 16,
    'chr18': 17,
    'chr19': 18,
    'chr20': 19,
    'chr21': 20,
    'chr22': 21,
    'chrX': 22,
    'chrY': 23,
    'chrM': 24
}

class Exon:
    def __init__(self, start, stop):
        self.start = start
        self.stop = stop

class Transcript:
    def __init__(self, id, type, start, stop):
        self.id = id
        self.type = type
        self.start = start
        self.stop = stop
        self.exons = []

class Gene:
    def __init__(self, id, name, type, chr, build, start, stop, strand):
        self.id = id
        self.name = name
        self.type = type
        self.chr = chr
        self.build = build
        self.start = start
        self.stop = stop
        self.strand = strand
        self.transcripts = {}

import sys

def reduce_gff3():
    genes = {}
    #read = 0
    for line in sys.stdin:
        if line.startswith('#'):
            continue
        #if (read%1000) == 0:
            #sys.stderr.write('Read {} lines\n'.format(read))
        #read += 1
        (chr, build, feat, start, stop, _, strand, _, attr) = line.strip().split('\t')
        if feat == 'gene':
            attrs = attr.split(';')
            gid, gname, gtype = "", "", ""
            # Find attribute section which matches 'gene_id'
            for attr in attrs:
                if attr.strip().startswith('gene_id'):
                    gid = attr.split()[1].replace('"', '')
                if attr.strip().startswith('gene_name'):
                    gname = attr.split()[1].replace('"', '')
                if attr.strip().startswith('gene_biotype'):
                    gtype = attr.split()[1].replace('"', '')
            #sys.stderr.write('Adding gene: {}({})\n'.format(gid, gname))
            genes[gid] = Gene(
                id=gid,
                name=gname,
                type=gtype,
                chr=chr,
                build=build,
                start=int(start),
                stop=int(stop),
                strand=strand
            )
        elif feat == 'transcript':
            attrs = attr.split(';')
            gid, txid, txtype = "", "", ""
            # Find attribute section
            for attr in attrs:
                if attr.strip().startswith('gene_id'):
                    gid = attr.split()[1].replace('"', '')
                elif attr.strip().startswith('transcript_id'):
                    txid = attr.split()[1].replace('"', '')
                elif attr.strip().startswith('transcript_biotype'):
                    txtype = attr.split()[1].replace('"', '')
            #sys.stderr.write('Adding tx: {}({}) - {}\n'.format(gid, gname, txid))
            genes[gid].transcripts[txid] = Transcript(
                id=txid,
                type=txtype,
                start=int(start),
                stop=int(stop)
            )
        elif feat == 'exon':
            attrs = attr.split(';')
            gid, txid = "", ""
            # Find attribute section
            for attr in attrs:
                if attr.strip().startswith('gene_id'):
                    gid = attr.split()[1].replace('"', '')
                elif attr.strip().startswith('transcript_id'):
                    txid = attr.split()[1].replace('"', '')
            #sys.stderr.write('    Adding exon: {}({}) - {}\n'.format(gid, gname, txid))
            genes[gid].transcripts[txid].exons.append(
                Exon(
                    start=int(start),
                    stop=int(stop)
                )
            )
        #else:
            #sys.stderr.write('Skipping {}\n'.format(feat))
    for gene in sorted(genes.values(), key=lambda g: (chrOrder[g.chr], g.start)):
        print('\t'.join([gene.chr, gene.build, 'gene', str(gene.start), str(gene.stop), '.', gene.strand, '.', 'gene_id "{}"; gene_name "{}"; gene_biotype "{}"'.format(gene.id, gene.name, gene.type)]))
        for tx in sorted(gene.transcripts.values(), key=lambda t: t.start):
            print('\t'.join([gene.chr, gene.build, 'transcript', str(tx.start), str(tx.stop), '.', gene.strand, '.', 'gene_id "{}"; gene_name "{}"; gene_biotype "{}"; transcript_id "{}"; transcript_biotype "{}"'.format(gene.id, gene.name, gene.type, tx.id, tx.type)]))
            for exon in sorted(tx.exons, key=lambda e: e.start):
                print('\t'.join([gene.chr, gene.build, 'exon', str(exon.start), str(exon.stop), '.', gene.strand, '.', 'gene_id "{}"; gene_name "{}"; gene_biotype "{}"; transcript_id "{}"; transcript_biotype "{}"'.format(gene.id, gene.name, gene.type, tx.id, tx.type)]))

if __name__ == '__main__':
    reduce_gff3()