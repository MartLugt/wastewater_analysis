#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="Calculate average mutation rate for fasta files")
    parser.add_argument('fastas', type=str, nargs='+', help='fasta files')
    parser.add_argument('-r', dest='ref', type=str, required=True, help='reference fasta')
    parser.add_argument('-o', dest='out_file', type=str, required=True, help='out file')
    args = parser.parse_args()

    match = []
    sub = []
    de = []
    ins = []

    out_dir = "/".join(args.out_file.split("/")[:-1])

    for filename in args.fastas:
        # Align with bbmap
        subprocess.check_call("bbmap.sh in={0} ref={1} out={2}/tmp.sam statsfile={2}/tmp.tsv"
            .format(filename, args.ref, out_dir), shell=True)
        
        # Read tsv with pandas
        df = pd.read_csv(out_dir+"/tmp.tsv", sep='\t', header=None, skiprows=17)

        print(df)

        match.append((float(df[3][0].strip().replace('%', ''))))
        sub.append((float(df[3][2].strip().replace('%', '')), int(df[4][2])))
        de.append((float(df[3][3].strip().replace('%', '')), int(df[4][3])))
        ins.append((float(df[3][4].strip().replace('%', '')), int(df[4][4])))

    sub_pct, sub_ct = zip(*sub)
    de_pct, de_ct = zip(*de)
    ins_pct, ins_ct = zip(*ins)    

    match_av = sum(match) / len(match)
    sub_av = (sum(sub_pct) / len(sub_pct), sum(sub_ct) / len(sub_ct))
    de_av = (sum(de_pct) / len(de_pct), sum(de_ct) / len(de_ct))
    ins_av = (sum(ins_pct) / len(ins_pct), sum(ins_ct) / len(ins_ct))

    print('match:', match_av)
    print("sub:", sub_av)
    print("de:", de_av)
    print("ins:", ins_av)

    file = open(args.out_file, "w")
    file.write('Error Type\tpct\tnum\n')
    file.write('Match:\t{}%\n'.format(match_av))
    file.write('Substitution:\t {}% \t{}\n'.format(sub_av[0], sub_av[1]))
    file.write('Deletion:\t {}% \t{}\n'.format(de_av[0], de_av[1]))
    file.write('Insertion:\t {}% \t{}\n'.format(ins_av[0], ins_av[1]))
    file.close()

    os.remove(out_dir+'/tmp.sam')
    os.remove(out_dir+'/tmp.tsv')


if __name__ == "__main__":
    sys.exit(main())