#!/usr/bin/env python3
"""
Module Docstring
"""

__author__ = "Harald Grove"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import time
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA as sklearnPCA
from adjustText import adjust_text
from scipy import stats

def transform(args):
    """ Transposes and cleans input matrix 
    """
    outfile = '{}.t.tsv'.format(args.transform.rsplit('.',1)[0])
    n = 2
    with open(args.transform, 'r') as fin:
        db = []
        for line in fin:
            l = line.strip().split()
            db.append(l[n:])
    names = []
    with open(outfile, "w") as fout:
        for line in zip(*db):
            fout.write("{}\n".format("\t".join(line)))

def transpose(args):
    """ Transposes input matrix 
    """
    outfile = '{}.t.tsv'.format(args.matrix.rsplit('.',1)[0])
    n = 1
    with open(args.matrix, 'r') as fin:
        db = []
        for line in fin:
            l = line.strip().split()
            db.append(l[n:])
    names = []
    with open(outfile, "w") as fout:
        for line in zip(*db):
            fout.write("{}\n".format("\t".join(line)))
            
def create_groups(args):
    if args.transform:
        infile = '{}.t.tsv'.format(args.transform.rsplit('.',1)[0])
    elif args.matrix is not None:
        infile = '{}.t.tsv'.format(args.matrix.rsplit('.',1)[0])
    else:
        infile = args.infile
    outfile = '{}.groups.tsv'.format(infile.rsplit('.', 1)[0])
    with open(infile, 'r') as fin, open(outfile, 'w') as fout:
        for line in fin:
            l = line.strip().split()
            fout.write('{}\t{}\n'.format(l[0], 'A'))

def dopca(args):
    if args.transform is not None:
        infile = '{}.t.tsv'.format(args.transform.rsplit('.',1)[0])
    elif args.matrix is not None:
        infile = '{}.t.tsv'.format(args.matrix.rsplit('.',1)[0])
    else:
        infile = args.infile
    if args.groups is not None:
        groupfile = args.groups
    else:
        groupfile = '{}.groups.tsv'.format(infile.rsplit('.', 1)[0])
    if args.header:
        df = pd.read_table(infile, header=0)
    else:
        df = pd.read_table(infile, header=None)
    extra = ""
    t = df.shape
    X = df.iloc[:, 1:].values
    df1 = pd.read_table(groupfile, header=None)
    names = df1.iloc[:,0]
    y1 = df1.iloc[:,1]
    group1 = df1.iloc[:,1].unique()
    ## Scale the dataset to unit scale (mean=0, variance=1).
    if args.standardize:
        from sklearn.preprocessing import StandardScaler
        X_std = StandardScaler().fit_transform(X)
        X = X_std
        extra = "std."
    comp = min(9, min(t)-1)
    sklearn_pca = sklearnPCA(n_components=comp)
    Y_sklearn = sklearn_pca.fit_transform(X)
    colorbar = ["#e6194b","#3cb44b","#ffe119","#0082c8",
            "#f58231","#911eb4","#46f0f0","#f032e6",
            "#d2f53c","#fabebe","#008080","#e6beff",
            "#aa6e28","#fffac8","#800000","#aaffc3",
            "#808000","#ffd8b1","#000080","#808080",
    ]
    colorbar2 = ["#FFFFFF", "#000000", "#aaaaaa"]
    with plt.style.context("seaborn-whitegrid"):
        texts = []
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        for lab, col in zip(group1, colorbar[0 : len(group1)]):
            plt.scatter(
                Y_sklearn[y1 == lab, 0],
                Y_sklearn[y1 == lab, 1],
                label=lab,
                c=col,
                s=100,
            )
        if args.labels:
            for ind, name in enumerate(names):
                #print(Y_sklearn[ind, 0], Y_sklearn[ind, 1], name)
                plt.text(Y_sklearn[ind, 0], Y_sklearn[ind, 1], name)
                #texts.append(plt.text(Y_sklearn[ind, 0], Y_sklearn[ind, 1], name))
        xp1, xp2 = sklearn_pca.explained_variance_[0:2]
        plt.xlabel(
            "PC 1 [{:.1f}%]".format(100 * xp1 / sum(sklearn_pca.explained_variance_))
        )
        plt.ylabel(
            "PC 2 [{:.1f}%]".format(100 * xp2 / sum(sklearn_pca.explained_variance_))
        )
        # use facecolors='none' for open circles
        # Add secondary grouping to the plot
        if False:
            for lab, col in zip(group2, colorbar2[0 : len(group2)]):
                plt.scatter(
                    Y_sklearn[y2 == lab, 0],
                    Y_sklearn[y2 == lab, 1],
                    s=10,
                    marker="o",
                    facecolors=col,
                    label=lab,
                    edgecolors=col,
                )
        if "NA" not in y1:
            plt.legend(loc=0, frameon=True)
        plt.tight_layout()
        #adjust_text(texts, arrowprops=dict(arrowstyle="->", color="red"))
        # plt.show()
        fig.savefig("{}.{}png".format(infile, extra))
    if not args.loadings:
        return
    # PCA loadings
    x = sklearn_pca.components_[0, :]
    y = sklearn_pca.components_[1, :]
    fig2, ax2 = plt.subplots(1, 1, figsize=(10, 10))
    plt.scatter(x, y)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    fig2.savefig("{}.loadings.{}png".format(args.infile,extra))
    #df2["pc1"] = x
    #df2["pc2"] = y
    #df2.to_csv(
    #    "{}_loadings.tsv".format(args.infile), header=True, index=False, sep="\t"
    #)


def main(args):
    """ Main entry point of the app """
    # If input file is the output from count2tpm, transform it first
    if args.transform is not None:
        transform(args)
    elif args.matrix is not None:
        transpose(args)
    # If no group file is specified, create one
    if args.groups is None:
        create_groups(args)
    dopca(args)
    if args.log:
        with open("README.txt", "a") as fout:
            fout.write("[{}]\t[{}]\n".format(time.asctime(), " ".join(sys.argv)))


if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    # Optional arguments with input
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
            "-i", "--infile",
            help="Pre-processed input file"
            )
    group.add_argument(
            "-t", "--transform",
            help="Count2tpm output file"
            )
    group.add_argument(
            "-m", "--matrix",
            help="Transpose matrix"
            )
    parser.add_argument(
            "-g", "--groups",
            help="Columns: sample name, group (no header)"
            )


    # Optional argument flag which defaults to False
    parser.add_argument(
            "-l","--log",
            action="store_false",
            default=True,
            help="Save command to 'README.txt'"
            )
    parser.add_argument(
            "-e", "--header",
            action="store_true",
            default=False,
            help="First row contains column headers"
            )
    parser.add_argument(
            "-b", "--labels",
            action="store_true",
            default=False,
            help="Add sample names to the plot"
            )
    parser.add_argument(
            "-a", "--loadings",
            action="store_true",
            help="Plot loadings"
            )
    parser.add_argument(
            "-s", "--standardize",
            action="store_true",
            help="Standardize the values"
            )

    # Optional verbosity counter (eg. -v, -vv, -vvv, etc.)
    parser.add_argument(
        "-v", "--verbose",
        action="count",
        default=0,
        help="Verbosity (-v, -vv, etc)"
    )

    # Specify output of '--version'
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__),
    )

    args = parser.parse_args()
    main(args)
