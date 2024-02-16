"""
Author: Ivy Strope
Created: 8/18/23
Contact: u247529@bcm.edu
Last Edited: 8/18/23 
"""
"""
Author: Ivy Strope
Created: 8/18/23
Contact: u247529@bcm.edu
Last Edited: 8/18/23
"""
#!/usr/bin/env python
import h5py
import argparse
import pandas as pd
import numpy as np
from scipy.sparse import csr_array

#SET UP HIERARCHICAL DATA STRUCTURE THAT MIMIC 10X OUTPUT
def create_hdf5(path,counts,organism,assembly):
    #Process Counts File
    data = pd.read_csv(counts,sep='\t',index_col='GENE')
    mtx = csr_array(data,dtype=int)  #convert to sparse format
    gene_names = list(data.index) #extract gene names
    ensembl = map_to_ensembl(genes = gene_names,organism = organism) #convert to ensembl id
    barcodes = list[data.columns.values] #extract spot barcode info
    feature_type = ['Gene_Expression'] * len(gene_names)

    #Create HDF5 File
    f = h5py.File(path,'w')
    #Create Hierarchical organization
    matrix = f.create_group('matrix')
    features = f.create_group('matrix/features')

    #Fill out matrix group
    matrix.create_dataset(name='barcodes',data=np.array(barcodes,dtype='S'))
    matrix.create_dataset(name='matrix',data=mtx.data,dtype='i')
    matrix.create_dataset(name='indices',data=mtx.indices,dtype='i')
    matrix.create_dataset(name='indptr',data=mtx.indptr,dtype='i')
    matrix.create_dataset(name='shape',data=mtx.shape,dtype='i')

    #Fill out features group
    features.create_dataset(name='_all_tag_keys',data = np.array('genome',dtype='S'))
    features.create_dataset(name = 'feature_type',
                                    data = np.array(feature_type,dtype='S'))
    features.create_dataset(name = 'genome',data = np.array(assembly,dtype='S'))
    features.create_dataset(name = 'id', data = np.array(ensembl,dtype='S'))
    features.create_dataset(name = 'name', data = np.array(gene_names,dtype='S'))

    f.close()


#MAP BETWEEN GENE SYMBOL AND ENSEMBL ID (Required for 10x)
def map_to_ensembl(genes,organism):

    if (organism == 'human'):
        gtf = pd.read_csv('species/human/ids.txt',sep='\t')
    else:
        gtf = pd.read_csv('species/mouse/ids.txt',sep='\t')
    gtf = gtf.reindex(genes)
    id = list(gtf['Geneid'])
    return id



parser = argparse.ArgumentParser(
        allow_abbrev=False,
        add_help=False)

parser.add_argument(
        '--counts',
        help='digital gene expression output file',
        type=str,
        required=True)

parser.add_argument(
        '--species',
        help='specify if mouse or human species',
        type=str,
        required=True)

parser.add_argument(
        '--assembly',
        help='specify assembly version',
        type=str,
        required=True)
parser.add_argument(
        '--output',
        help='output file path',
        type=str,
        required=True)


#RUN SCRIPTS
args = parser.parse_args()
create_hdf5(path = args.output,
        counts = args.counts,
        organism = args.species,
        assembly = args.assembly)
