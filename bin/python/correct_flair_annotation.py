from sys import stdin

def read_ann(path):
    io = open(path, 'r') if path != "-" else stdin
    with io as f:
        last_gene = ["", 0, 0, "", ""]
        transcripts = []
        dat = f.readline().strip().split('\t')
        while not dat[2] == "transcript":
            dat = f.readline().strip().split('\t')
        gene_id = [attr.split('"')[1].strip() for attr in dat[8].split(';') if attr.strip().startswith("gene_id")][0]
        last_gene = [dat[0], int(dat[3]), int(dat[4]), dat[6], gene_id]
        transcripts = ['\t'.join(dat)]
        for line in f:
            dat = line.strip().split('\t')
            if (dat[2]=="transcript") & (last_gene[4]!=gene_id):
                yield ["{}\tFLAIR\tgene\t{}\t{}\t.\t{}\t.\tgene_id \"{}\";".format(last_gene[0], last_gene[1], last_gene[2], last_gene[3], last_gene[4])]+transcripts
                last_gene = [dat[0], int(dat[3]), int(dat[4]), dat[6], gene_id]
                transcripts = [line.strip()]
            else:
                transcripts.append(line.strip())
        yield ["{}\tFLAIR\tgene\t{}\t{}\t.\t{}\t.\tgene_id \"{}\";".format(last_gene[0], last_gene[1], last_gene[2], last_gene[3], last_gene[4])]+transcripts

def main():
    from sys import argv
    annotation_path = argv[1]
    for geneset in read_ann(annotation_path):
        for ln in geneset:
            print(ln)
main()