import scanpy as sc
from sklearn.cluster import KMeans
import argparse

def arguments():
    parser = argparse.ArgumentParser(description='Process Anndata Objects')
    
    parser.add_argument('--input', help='input anndata object')
    parser.add_argument('--output', help='output processed anndata object')
    return parser.parse_args()

def Downstream(input,output):
    adata = sc.read_h5ad(input)
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith(("MT-","mt-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    #Detect Tissue Segments
    metrics = adata.obs[['n_genes_by_counts','log1p_n_genes_by_counts','total_counts','log1p_total_counts']]
    kmeans = KMeans(2,n_init = 1000)
    kmeans.fit(metrics)
    adata.obs['in_tissue_pred'] = kmeans.labels_

    #Normalize
    sc.pp.normalize_total(adata,inplace=True)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata,inplace=True)
    #principal component analysis,clustering,low dimensional projections
    sc.tl.pca(adata)
    sc.pp.neighbors(adata,n_neighbors=10,n_pcs=40)
    sc.tl.umap(adata)

    sc.tl.leiden(adata)
    adata.write_h5ad(output)


args = arguments()
Downstream(args.input,args.output)
