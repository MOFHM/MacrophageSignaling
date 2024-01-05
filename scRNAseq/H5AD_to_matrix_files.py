os.chdir(dir)
import scanpy as sc
from scipy import io
adata = sc.read_h5ad("file.h5ad")
adata = adata.raw.to_adata()

with open('barcodes.tsv', 'w') as f:
    for item in adata.obs_names:
        f.write(item + '\n')

with open('features.tsv', 'w') as f:
    for item in ['\t'.join([x,x,'Gene Expression']) for x in adata.var_names]:
        f.write(item + '\n')

io.mmwrite('matrix', adata.X.T)

adata.obs.to_csv('metadata.csv')
