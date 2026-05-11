import anndata
import pandas as pd
from ALLCools.mcds import MCDS

mcds_path = "genebody_frac.mcds/"
adata_path = "adata.with_coords_final_pooled.h5ad"

adata = anndata.read_h5ad(adata_path)
adata = adata[adata.obs["Group"].isin(["Control", "Stressed"])].copy()

gene_mcds = MCDS.open(mcds_path, use_obs=adata.obs_names)

da = gene_mcds["genebody_da_frac"]  # (cell, genebody, mc_type)

for mc_type in ["CGN", "CHN"]:
    mat = da.sel(mc_type=mc_type).transpose("genebody", "cell").to_pandas()
    mat.to_csv(f"genebody_frac_{mc_type}.csv")
    print("Wrote:", f"genebody_frac_{mc_type}.csv", mat.shape)

meta = adata.obs[["Group"]].copy()
meta["sample"] = meta.index
meta = meta[["sample", "Group"]]
meta.to_csv("sample_metadata.csv", index=False)
print("Wrote: sample_metadata.csv", meta.shape)
