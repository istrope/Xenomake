import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
import argparse

def arguments(required=True):
    parser=argparse.ArgumentParser(
        allow_abbrev=False,
        add_help=False
    )
    parser.add_argument(
        '--counts',
        help='digital gene expression output file',
        type=str,
        required=True
    )
    parser.add_argument(
        '--spatial_coordinates',
        help='path to spatial barcode file',
        type=str,
        required=True
    )
    parser.add_argument(
        '--output',
        help='output file path',
        type=str,
        required=True)

args = arguments()
expr = pd.read_csv(args.counts,sep='\t',index_col='GENE')
expr = expr.transpose()
cell_ids = expr.index

#KMeans Cluster of 2 (predict in and out of tissue)
kmeans = KMeans(2)
kmeans.fit(expr)
in_tissue=kmeans.fit_predict(expr)
in_tissue=pd.DataFrame(in_tissue,index=cell_ids,columns='in_tissue')

#import barcode spatial information
spatial_barcode=pd.read_csv(args.spatial_coordinates)
spatial_barcode.columns = ['barcode','x_pos','y_pos']
spatial_barcode.index = spatial_barcode['barcode']

#Format Tissue Positions List
tissue_file = spatial_barcode.merge(in_tissue)
tissue_file = tissue_file[['barcode','in_tissue','x_pos','y_pos']]
tissue_file.to_csv(args.output,header=None)

