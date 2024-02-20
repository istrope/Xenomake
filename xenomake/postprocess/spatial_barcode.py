import pandas as pd
from anndata import AnnData
import scanpy as sc
import numpy as np
from sklearn.cluster import KMeans
import argparse
from sklearn.neighbors import KDTree

#KMeans Cluster of 2 (predict in and out of tissue)
def pred_tissue(expr):
    norm_expr = sc.pp.normalize_total(AnnData(np.array(expr)),target_sum=1,inplace=False)['X']
    log_expr = sc.pp.log1p(norm_expr)
    #First prediction of in-tissue based upon umi counts
    kmeans = KMeans(n_clusters=2,n_init = 10)
    kmeans.fit(log_expr)
    in_tissue=kmeans.labels_
    
    return in_tissue

def arguments():
    parser=argparse.ArgumentParser(allow_abbrev=False,add_help=False)
    parser.add_argument('--counts',help='digital gene expression output file',type=str,required=True)
    parser.add_argument('--spatial_coordinates',help='path to spatial barcode file',type=str,required=True)
    parser.add_argument('--output',help='output file path',type=str,required=True)

    return parser.parse_args()


def read_data(counts,coordinates):
    expr = pd.read_csv(counts,sep='\t',index_col='GENE')
    expr = expr.transpose()
    spatial_barcode=pd.read_csv(coordinates,index_col='cell_bc')
    return expr, spatial_barcode



args = arguments()

expr,spatial_barcode = read_data(args.counts,args.spatial_coordinates)
#Create Dataframe for subsequent merging
in_tissue=pd.DataFrame(list(pred_tissue(expr)),expr.index)
in_tissue.columns = ['KMeans']

#Format Tissue Positions List
tissue_file = spatial_barcode.merge(in_tissue,left_index=True,right_index=True)
tissue_file = tissue_file[['KMeans','x_pos','y_pos']]

#Write out file
tissue_file.to_csv(args.output)
