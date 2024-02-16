import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sb
import argparse

def arguments():
    parser=argparse.ArgumentParser(
        allow_abbrev=False,
        add_help=False)
    
    parser.add_argument(
                '--input',
                help='name of input in h5ad format',
                type=str,
                required=True)
    parser.add_argument(
                '--out',
                help='name of output png file',
                type=str,
                required=True)
    return parser.parse_args()


def qc(h5ad,out):
    adata = sc.read_h5ad(h5ad)
    adata.var['mt'] = adata.var_names.str.startswith(('MT-','mt-'))
    sc.pp.calculate_qc_metrics(adata,inplace=True,qc_vars=['mt'])

    fig, axs = plt.subplots(1,4, figsize=(15,4))
    fig.suptitle('Covariates for filtering')
    sb.distplot(adata.obs['total_counts'], kde=False, ax = axs[0])
    sb.distplot(adata.obs['total_counts'][adata.obs['total_counts']<10000], kde=False, bins=40, ax = axs[1])
    sb.distplot(adata.obs['n_genes_by_counts'], kde=False, bins=60, ax = axs[2])
    sb.distplot(adata.obs['n_genes_by_counts'][adata.obs['n_genes_by_counts']<4000], kde=False, bins=60, ax = axs[3])

    count_out = out + '_qc.png'
    fig.savefig(count_out)
    plt.close(fig)

args = arguments()
qc(h5ad=args.input,out=args.out)

