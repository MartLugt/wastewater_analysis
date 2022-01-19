#!/usr/bin/env python3

import sys
import argparse
import subprocess
import pandas as pd
import json
from math import log10, floor

def main():
    parser = argparse.ArgumentParser(description="Analyse kallisto output")
    parser.add_argument('predictions', type=str, nargs='+')
    parser.add_argument('--metadata', dest='metadata', type=str, required=True, help="metadata file")
    parser.add_argument('--output', dest='output', type=str, required=True, help="output file")
    parser.add_argument('--voc', dest='voc', type=str, help="comma-separated list of strains of interest, output abundance for these only")
    args = parser.parse_args()

    voc_list = args.voc.split(',')

    averages = dict()

    for voc in voc_list:
        averages[voc] = (0, 0, 0, 0, 0)

    for pred in args.predictions:
        voc = pred.split('/')[-2].split('_')[0]

        outfilepath = pred.split('/')
        outfilepath[-1] = "predictions_other.tsv"
        outfile = "/".join(outfilepath)

        subprocess.check_call("python pipeline/output_abundances.py --metadata {0} --voc {3} -o {1} {2}"
            .format(args.metadata, outfile, pred, args.voc), shell=True)

        df = pd.read_csv(outfile, sep='\t', header=None, skiprows=3)    

        print(df)

        other = 0
        B117 = 0
        own = 0

        for row in range(len(df.index)):
            if (df[0][row] == voc):
                own = df[2][row]
                continue
            if (df[0][row] == "B.1.1.7"):
                B117 = df[2][row]
                continue
            if (df[0][row] == "other"):
                other = df[2][row]
                continue

        run_info = pred.split('/')
        run_info[-1] = "run_info.json"
        run_info = "/".join(run_info)

        f = open(run_info)

        data = json.load(f)
        pseudo = data["p_pseudoaligned"]

        f.close()

        av = averages[voc]
        averages[voc] = (av[0] + 1, av[1] + pseudo, av[2] + own, av[3] + B117, av[4] + other)

    f = open(args.output, "w")
    f.write('Variant\tpct pseudoaligned\tpct correct\tpct b.1.1.7\tpct other\n')

    for voc in voc_list:
        avc = averages[voc]
        av = (avc[1] / avc[0], avc[2] / avc[0], avc[3] / avc[0], avc[4] / avc[0])
        print(voc, av)

        f.write('{}\t{}\t{}\t{}\t{}\n'.format(voc, round_sig(av[0]), round_sig(av[1]), round_sig(av[2]), round_sig(av[3])))

    f.close()


def round_sig(x, sig=3):
    if x == 0:
        return x
    return round(x, sig - int(floor(log10(abs(x)))) - 1)


if __name__ == "__main__":
    sys.exit(main())

