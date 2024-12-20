def main():
    from sys import argv,stdin
    dcs_seqs = []
    with open(argv[0]) as dcs:
        for line in dcs:
            dcs_seqs.append(line.strip())
    with stdin as fastq:
        for line in fastq:
            if (line[0] == '@') and (line[1:] not in dcs_seqs):
                print(line.strip())
                for _ in [1,2,3]: print(fastq.readline().strip())

main()