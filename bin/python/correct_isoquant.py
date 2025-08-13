def read_ann(path):
    with open(path, 'r') as f:
        for line in f:
            if 'novel' in line:
                gene_id = [attr.split('"')[1].strip() for attr in line.strip().split('\t')[-1].split(';') if attr.strip().startswith("gene_id")][0]
                line = line.strip() + " gene_name \"{}\"; gene_type \"novel\";".format(gene_id)
            yield line.strip()

def main():
    from sys import argv
    annotation_path = argv[1]
    for line in read_ann(annotation_path):
        print(line)
main()