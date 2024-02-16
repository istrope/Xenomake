import scanpy as sc
import pandas as pd
import argparse

parser=argparse.ArgumentParser(
        allow_abbrev=False,
        add_help=False)

parser.add_argument(
                '--counts',
                help='digital gene expression output file',
                type=str,
                required=True)

parser.add_argument(
                '--spatial_coordinates',
                help='path to spatial barcode file',
                type=str,
                required=True)

parser.add_argument(
                '--output',
                help='output file path',
                type=str,
                required=True)




args = parser.parse_args()


adata = sc.read(args.counts)
adata = adata.transpose()
adata.var_names_make_unique()
positions = pd.read_csv(args.spatial_coordinates,header=None,index_col = 0)
positions.columns = [
              'barcode',
              'x_pos',
              'y_pos']
adata.obs = adata.obs.merge(positions, how='left',left_index=True,right_index=True)
adata.obsm['spatial'] = adata.obs[['x_pos', 'y_pos']].to_numpy()

adata.write_h5ad(args.output)
