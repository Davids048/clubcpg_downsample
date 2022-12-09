import argparse
import numpy as np
import pandas as pd
import os

def get_cluster_patterns(clubcpgPathLst: list):
    allPattern = pd.DataFrame(columns=["cpg_number","cpg_pattern", "methylation"])
    for path in clubcpgPathLst:
        if path.endswith(".csv"):
            clubcpg = pd.read_csv(path)
            clubcpg = clubcpg[["cpg_number","cpg_pattern", "methylation"]]
            clubcpg = clubcpg.drop_duplicates(inplace=False).reset_index(drop = True)
            allPattern = pd.concat([allPattern,clubcpg], axis=0)
            allPattern = allPattern.drop_duplicates(inplace = False).reset_index(drop = True)
        else:
            print(path + "\ not a csv file")
            continue

    return allPattern


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # parser.add_argument("-A", help="absolute path to the first clubcpg cluster output", default=None)
    # parser.add_argument("-B", help="absolute path to the second clubcpg cluster output", default=None)
    parser.add_argument("--samples", help="absolute path to the  clubcpg cluster outputs", default=None, nargs="+")
    parser.add_argument("-o", "--output", help="folder to save imputed coverage data", default=None)
    parser.add_argument("-name", "--name", help =" file name to save to, start with \"/\"", default = "/lc_reads1.csv")
    args = parser.parse_args()

    if not args.output:
        output_folder = os.path.dirname(args.A + "/../DEC7")
    else:
        output_folder = args.output
    try:
        os.mkdir(output_folder)
    except FileExistsError:
        print("Output folder already exists... no need to create it...")

    patterns = get_cluster_patterns(list(args.samples))
    patterns.to_csv(args.output + args.name)
    print(patterns)
    print("done")