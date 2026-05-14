def attrs_to_dict(attrs):
    '''Extract attributes from GTF/GFF3 annotation.'''
    attr_dict = {}
    for attr in attrs.split(';'):
        if attr.strip() == '':
            continue
        key, value = attr.strip().split(' ')
        attr_dict[key] = value.strip('"')
    return attr_dict

def read_ann(path):
    '''Add gene_name + gene_type to IsoQuant GTF/GFF3 transcript annotations.'''
    with open(path, 'r') as f:
        curr_gene_info = {}
        for ann in f:
            if ann.startswith('#'):
                yield ann.strip()
                continue
            if ann.strip().split('\t')[2] == 'gene':
                curr_gene_info = attrs_to_dict(ann.strip().split('\t')[-1])
                if (not curr_gene_info.get('gene_name')): 
                    gene_id = curr_gene_info.get('gene_id', '')
                    gene_name = gene_id.replace('novel_gene_', 'Novel gene ')
                    curr_gene_info['gene_name'] = gene_name
                    curr_gene_info['gene_biotype'] = 'Novel gene'
                    yield ann.strip() + ' gene_name \"{}\"; gene_biotype \"Novel gene\"'.format(gene_name)
                yield ann.strip()
                continue
            if ann.strip().split('\t')[2] == 'transcript':
                yield ann.strip() + ' gene_name \"{}\"; gene_biotype \"{}\"'.format(curr_gene_info.get('gene_name', ''), curr_gene_info.get('gene_biotype', ''))
                continue
            yield ann.strip()

def main():
    from sys import argv
    annotation_path = argv[1]
    for line in read_ann(annotation_path):
        print(line)
main()
