def main():
    from sys import stdin
    novel_genes = set()
    for line in stdin:
        dat = line.strip().split('\t')
        if dat[2] == "transcript":
            gene_id = [attr.split('"')[1].strip() for attr in dat[8].split(';') if attr.strip().startswith("gene_id")][0]
            transcript_id = [attr.split('"')[1].strip() for attr in dat[8].split(';') if attr.strip().startswith("transcript_id")][0]
            gene_name = "Novel_"+gene_id
            novel_genes.add(gene_name)
            print(gene_id,transcript_id,"Novel Isoform",gene_name,sep=',')
main()