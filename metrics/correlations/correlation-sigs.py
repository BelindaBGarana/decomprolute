#!/usr/local/bin/python

import pandas as pd
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--proteomics1', dest='profile1',
                        help='Deconvoluted matrix of proteomics data using signature 1')
    parser.add_argument('--proteomics2', dest='profile2',
                        help='Deconvoluted matrix of proteomics data using signature 2')
    parser.add_argument('--spearOrPears', dest='sop', help='Use spearman or pearson correlation',
                        default='pearson')
    # parser.add_argument('--output', dest='output',
    #                     help='Output file for the correlation values')
    opts = parser.parse_args()

    pro1 = pd.read_csv(opts.profile1, sep='\t', index_col=0)
    pro2 = pd.read_csv(opts.profile2, sep='\t', index_col=0)

    pro1Cols = list(pro1.columns)
    pro2Cols = list(pro2.columns)
    pro1Rows = list(pro1.index)
    pro2Rows = list(pro2.index)
    intersectCols = list(set(pro1Cols) & set(pro2Cols))
    intersectRows = list(set(pro1Rows) & set(pro2Rows))

    pro1 = pro1.loc[intersectRows, intersectCols]
    pro2 = pro2.loc[intersectRows, intersectCols]

    if opts.sop == 'pearson':
        corrList = [pro1[sample].corr(pro2[sample]) for sample in intersectCols]
    else:
        corrList = [pro1[sample].corr(
            pro2[sample], method='spearman') for sample in intersectCols]
    correlations = pd.Series(corrList)
    correlations.index = intersectCols
    correlations.to_csv("corr.tsv", sep='\t', header=False)


if __name__ == '__main__':
    main()
