'''
Author: Ivy Strope
Date: 6/7/23
Contact: u247529@bcm.edu
'''
 #Determine how to change this based upon input
sample=config['sample']
OUTDIR=config['outdir']
spatial_dir=

import scanpy as sc
import json
import cv2
import os
import pandas as pd
import numpy as np

rule scanpy_load:
    input:
        count_matrix=config['count_matrix']
        spatial_path=config['spatial_dir']
    params:
        sample=config['sample']
    output:
        '{sample}_{organism}_anndata'
    run:
        adata = sc.read(count_matrix)
        adata.var_names_make_unique()
        positions = pd.read_csv(os.path.join(spatial_path,"tissue_positions_list.csv"),header=None)
        positions.columns = [
              'barcode',
              'in_tissue',
              'array_row',
              'array_col',
              'pxl_col_in_fullres',
              'pxl_row_in_fullres',
          ]
          positions.index = positions['barcode']
    #positions['pxl_row_in_fullres'],positions['pxl_col_in_fullres'] = transform(positions['pxl_col_in_fullres'],positions['pxl_row_in_fullres'],name)
        adata.obs = adata.obs.join(positions, how="left")
        adata.obsm['spatial'] = adata.obs[['pxl_row_in_fullres', 'pxl_col_in_fullres']].to_numpy()

    #Filter Genes
        sc.pp.filter_genes(adata,min_counts = 250)
        sc.pp.filter_genes(adata,min_cells=25)

        sc.pp.filter_cells(adata,)
    #Filter Cells by Necrosis and inTissue

        necrosis = pd.read_csv(necrosis_file)
        necrosis.columns = ['barcode']
        keep = positions[~positions.barcode.isin(necrosis.barcode)]
        keep = keep.loc[keep['in_tissue'] ==1]
        keep_cells = keep.barcode

        adata = adata[keep_cells]
    #adata.obs.drop(columns=['barcode', 'pxl_row_in_fullres', 'pxl_col_in_fullres'], inplace=True)

rule process_scanpy:
    input:
        '{sample}_{organism}_anndata'
    output:
        '{sample}_{organism}_anndata'
    run:
        #Normalization,log transform, hvgs
            sc.pp.normalize_total(adata, target_sum=1e4,inplace=True)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata,inplace=True)
        #principal component analysis,clustering,low dimensional projections
            sc.tl.pca(adata)
            sc.pp.neighbors(adata,n_neighbors=10,n_pcs=30)
            sc.tl.umap(adata)
            sc.tl.leiden(adata,resolution = 0.4)
            #save file
            adata.write({spacemake.output})

rule add_image:
    input:
    output:
    run:
        import cv2
        import os
        img = cv2.imread(os.path.join(spatial_dir,'tissue_hires_image.png'))
        import json
        with open(os.path.join(spatial_dir,'scalefactors_json.json')) as user_file:
            file_contents = user_file.read()
        scalefactors = json.loads(file_contents)

        spatial_key = "spatial"
        library_id = "tissue"

        adata.uns[spatial_key] = {library_id: {}}
        adata.uns[spatial_key][library_id]["images"] = {}
        adata.uns[spatial_key][library_id]["images"] = {"hires": img}
        adata.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": scalefactors['tissue_hires_scalef'],
                                                            'fiducial_diameter_fullres': scalefactors['fiducial_diameter_fullres'],
                                                            "spot_diameter_fullres": scalefactors['spot_diameter_fullres']}

rule commot:
    input:
        '{sample}_{organism}_anndata'
    output:
        '{sample}_{organism}_anndata'
    run:
        import commot as ct
        df_cellchat = ct.pp.ligand_receptor_database(species='human', signaling_type='Secreted Signaling', database='CellChat')
        df_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, adata, min_cell_pct=0.05)
        ct.tl.spatial_communication(adata, database_name='cellchat', df_ligrec=df_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True)
        adata.write_h5ad(filename={snakemake.output},compression='gzip')

rule DE_genes:
    input:
    output:
    run:

rule gsea:
    input:
    output:
    run:
