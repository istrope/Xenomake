import pandas as pd
import argparse

def arguments():
    parser=argparse.ArgumentParser(allow_abbrev=False,add_help=False)
    parser.add_argument('--counts',help='digital gene expression output file',type=str,required=True)
    parser.add_argument('--output',help='output file path',type=str,required=True)

args = arguments()

df = pd.read_csv(args.counts,sep='\t')
cols = list(df.columns[1:])
long =pd.melt(df,id_vars='GENE',value_vars=cols,var_name='cell_ID',value_name='Count')
drop_zeros=long.loc[long['Count'] !=0,]

drop_zeros.to_csv(args.output,index=False)