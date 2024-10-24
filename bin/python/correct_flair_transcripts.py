def main():
    from sys import argv
    from re import search
    transcript_path = argv[1]
    with open(transcript_path, 'r') as fa:
        for line in fa:
            if line.startswith('>'):
                if ('range=' in line): # repeat
                    lidx = search(r"=[\+|-]", line.strip()).end()
                    print(line.strip()[:(lidx+1)])
                else:
                    tx_id = '_'.join(line.strip().split('_')[:-1])
                    print(tx_id)
            else:
                print(line.strip())

main()