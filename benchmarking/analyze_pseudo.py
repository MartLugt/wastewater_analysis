#!/usr/bin/env/python3

import sys
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Analyse misprediction from sam file")
    parser.add_argument('sam', type=str)
    parser.add_argument('--metadata', type=str, required=True)
    parser.add_argument('--voc', type=str)
    parser.add_argument('--output', type=str)
    args = parser.parse_args()

    vocs = args.voc.split(',')

    meta = pd.read_csv(args.metadata, sep='\t', header=0, dtype=str)

    # match = []

    count = 0
    res = dict()

    vocs.append("other")
    vocs.append("None")

    for voc in vocs:
        res[voc] = dict()
        for voc2 in vocs:
            res[voc][voc2] = 0

    with open(args.sam) as f:
        for line in f:
            if line[0] == '@':
                continue

            parts = line.split()

            read_virus = '-'.join(parts[0].split('-')[:-1]).split('|')[0]
            trans_virus = parts[2].split('|')[0]

            # print(read_virus, "->", trans_virus)

            lineage = meta.loc[meta["Virus name"] == read_virus]["Pango lineage"]
            read_pango = lineage.iloc[0]

            lineage = meta.loc[meta["Virus name"] == trans_virus]["Pango lineage"]
            trans_pango = lineage.iloc[0]

            # print(read_pango, "->", trans_pango)

            # match.append((read_pango, trans_pango))

            v0 = read_pango if read_pango in vocs else 'other'
            v1 = trans_pango if trans_pango in vocs else 'other'

            res[v0][v1] += 1

            count += 1
            if count % 1000 == 0:
                print("Reads read:", count, end="\r", flush=True)

    # for m in match:
    #     v0 = m[0] if m[0] in vocs else 'other'
    #     v1 = m[1] if m[1] in vocs else 'other'

    #     res[v0][v1] += 1
    print("\n")

    with open(args.output, 'w') as f:
        header = "\t" + "\t".join(vocs)

        print(header)
        header += '\n'
        f.write(header)

        for voc in vocs:
            line = voc + "\t"

            for voc2 in vocs:
                line += str(res[voc][voc2]) + '\t'

            print(line)
            line += '\n'
            f.write(line)



if __name__ == "__main__":
    sys.exit(main())