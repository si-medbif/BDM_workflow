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


def dopca(args):
    if args.header:
        df = pd.read_table(args.infile, header=0)
    else:
        df = pd.read_table(args.infile, header=None)
    if args.groups is not None:
        df1 = pd.read_table(args.groups, header=0)
        category1 = "NA"
        y1 = np.array(["NA"] * (t[0]))
        group1 = ["NA"]
        group2 = ["NA"]
    else:
        t = df.shape
        X = df.iloc[:, 1:].values
        category1 = "NA"
        y1 = np.array(["NA"] * (t[0]))
        group1 = ["NA"]
        group2 = ["NA"]
    # y1 is a list with sample names, or NA.
    if args.varlabels is not None:
        df2 = pd.read_table(args.varlabels, header=0)
    ## Scale the dataset to unit scale (mean=0, variance=1).
    # from sklearn.preprocessing import StandardScaler
    # X_std = StandardScaler().fit_transform(X)
    comp = min(9, min(t)-1)
    sklearn_pca = sklearnPCA(n_components=comp)
    Y_sklearn = sklearn_pca.fit_transform(X)
    colorbar = [
        "#e6194b",
        "#3cb44b",
        "#ffe119",
        "#0082c8",
        "#f58231",
        "#911eb4",
        "#46f0f0",
        "#f032e6",
        "#d2f53c",
        "#fabebe",
        "#008080",
        "#e6beff",
        "#aa6e28",
        "#fffac8",
        "#800000",
        "#aaffc3",
        "#808000",
        "#ffd8b1",
        "#000080",
        "#808080",
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
            for ind, name in enumerate(y1):
                print(Y_sklearn[ind, 0], Y_sklearn[ind, 1], name)
                #plt.text(Y_sklearn[ind, 0], Y_sklearn[ind, 1], name)
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
        fig.savefig("{}.png".format(args.infile))
    if args.varlabels is None:
        return
    # PCA loadings
    x = sklearn_pca.components_[0, :]
    y = sklearn_pca.components_[1, :]
    fig2, ax2 = plt.subplots(1, 1, figsize=(10, 10))
    plt.scatter(x, y)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    fig2.savefig("{}.loadings.png".format(args.infile))
    df2["pc1"] = x
    df2["pc2"] = y
    df2.to_csv(
        "{}_loadings.tsv".format(args.infile), header=True, index=False, sep="\t"
    )


def main(args):
    """ Main entry point of the app """
    dopca(args)
    if args.log:
        with open("README.txt", "a") as fout:
            fout.write("[{}]\t[{}]\n".format(time.asctime(), " ".join(sys.argv)))


if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    # Required positional argument
    parser.add_argument("infile", help="Input file")

    # Optional argument flag which defaults to False
    parser.add_argument(
        "-l",
        "--log",
        action="store_true",
        default=True,
        help="Save command to 'README.txt'",
    )
    parser.add_argument(
            "-e", "--header", action="store_true", default=False, 
            help="First row contains column headers"
    )
    parser.add_argument(
            "-n", "--samples", action="store_true", default=False, 
            help="First column contains sample names"
    )
    parser.add_argument(
            "-t", "--transpose", action="store_true", default=False,
            help="Transpose matrix"
    )
    parser.add_argument(
        "-b", "--labels", action="store_true", default=False, help="Add labels to plot"
    )

    # Optional argument which requires a parameter (eg. -d test)
    parser.add_argument("-g", "--groups", action="store", help="Sample information")
    parser.add_argument(
        "-a", "--varlabels", action="store", help="Variable information"
    )

    # Optional verbosity counter (eg. -v, -vv, -vvv, etc.)
    parser.add_argument(
        "-v", "--verbose", action="count", default=0, help="Verbosity (-v, -vv, etc)"
    )

    # Specify output of '--version'
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__),
    )

    args = parser.parse_args()
    main(args)
