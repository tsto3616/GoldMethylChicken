# Adpated from the poly-CpG package/ GitHub repository from the same Author of the GitHub repository
import numpy as np
import pandas as pd

# Define standard GTF column names
gtf_cols = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]

# Read the file
gtf = pd.read_csv(
    "D:/bird_methylome/CHK/galGal6_gtf.gz", 
    sep="\t", 
    comment="#", 
    header=None, 
    names=gtf_cols
) 

# below is to parse the gtf file through python - this line extracts only gene level data
genes_annot = gtf[gtf["feature"] == "CDS"].copy()

genes_annot["gene_symbol"] = genes_annot["attribute"].str.extract(r'transcript_id "([^"]+)"')

# 1. Group and calculate absolute physical boundaries
compressed_cds = (
    genes_annot.groupby(["seqname", "strand", "gene_symbol"])
    .agg(
        # Always find the absolute outer numerical coordinates first
        physical_start=("start", "min"),
        physical_end=("end", "max"),
    )
    .reset_index()
)

# 2. Add flexible logic for biological directionality
compressed_cds["biological_start"] = np.where(
    compressed_cds["strand"] == "+",
    compressed_cds["physical_start"],  # On + strand, biological start is the lowest
    compressed_cds["physical_end"],  # On - strand, biological start is the highest
)

compressed_cds["biological_end"] = np.where(
    compressed_cds["strand"] == "+",
    compressed_cds["physical_end"],  # On + strand, biological end is the highest
    compressed_cds["physical_start"],  # On - strand, biological end is the lowest
)

# just to check that the result is ok
print(compressed_cds["gene_symbol"].head(30))

compressed_cds["chromosome"] = (
    compressed_cds["seqname"]
    .astype(str)
    .str.replace('"', '', regex=False)
    .str.replace("seqname", "", regex=False)
    .str.strip()
)

# just to check that the result is ok
print(compressed_cds["chromosome"].head(5))

# TSS - TSS is the transcription start site
compressed_cds["tss"] = np.where(
    compressed_cds["strand"] == "+",
    compressed_cds["biological_start"],
    compressed_cds["biological_end"]
)

# basal regulatory domain - important in the GREAT method
compressed_cds["reg_start"] = compressed_cds["tss"] - 5000
compressed_cds["reg_end"]   = compressed_cds["tss"] + 1000

# sort by chr, tss
compressed_cds = compressed_cds.sort_values(["chromosome", "tss"])

# extend to midpoints between adjacent genes - so as to prevent the uneven assignment of CpGs to 
# particular genes 
for chr_name, df_chr in compressed_cds.groupby("chromosome"):
    idxs = df_chr.index.to_list()
    tss = df_chr["tss"].values
    if len(tss) > 1:
        mids = (tss[:-1] + tss[1:]) / 2
        compressed_cds.loc[idxs[0], "reg_end"] = mids[0]
        for i in range(1, len(df_chr)-1):
            compressed_cds.loc[idxs[i], "reg_start"] = mids[i-1]
            compressed_cds.loc[idxs[i], "reg_end"]   = mids[i]
        compressed_cds.loc[idxs[-1], "reg_start"] = mids[-1]

# load the file for the gene annotation
annot = pd.read_csv("D:/bird_methylome/CHK/chick_galGal6/chicken_galGal6_merged.csv")

matches = []

for (chromosome, strand), pos_grp in annot.groupby(["chromosome", "strand"]):

    ranges = compressed_cds[
        (compressed_cds["chromosome"] == chromosome)
        & (compressed_cds["strand"] == strand)
    ]

    if ranges.empty:
        continue

    for _, pos_row in pos_grp.iterrows():

        matched_ranges = ranges[
            (ranges["reg_start"] <= pos_row["position"])
            & (ranges["reg_end"] >= pos_row["position"])
        ]

        if matched_ranges.empty:
            continue

        merged = pd.concat(
            [pd.DataFrame([pos_row] * len(matched_ranges)).reset_index(drop=True),
             matched_ranges.reset_index(drop=True)],
            axis=1
        )

        matches.append(merged)

print(matches)

result = pd.concat(matches, ignore_index=True)

result.to_csv("D:/bird_methylome/CHK/chick_galGal6/chick_galGal6_merged_with_gene.csv", index=False)
