"""
This script is intended to be used with Cell Cell Interactions Downstream of analysis as described on the Github README.md
Conver_long.py takes a UMI count matrix and outputs a new matrix with three columns (Gene Name, Cell Barcode, UMI Count). It also removes all zero counts from the file
"""
import pandas as pd
import argparse

def arguments():
    parser=argparse.ArgumentParser(allow_abbrev=False,add_help=False)
    parser.add_argument('--counts',help='digital gene expression output file',type=str,required=True)
    parser.add_argument('--output',help='output file path',type=str,required=True)
    return parser.parse_args()

args = arguments()

df = pd.read_csv(args.counts,sep='\t')
cols = list(df.columns[1:])
long =pd.melt(df,id_vars='GENE',value_vars=cols,var_name='cell_ID',value_name='Count')
drop_zeros=long.loc[long['Count'] !=0,]

drop_zeros.to_csv(args.output,index=False)
