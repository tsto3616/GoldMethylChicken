import pandas as pd
from Bio import SeqIO
from scipy.stats import hypergeom

CHK_1_H = pd.read_csv("data/chicken_galGal6/CHK_1_H_merged_with_gene.csv")
IE_1_H = pd.read_csv("data/IE_galGal6/IE_1_H_merged_with_gene.csv")

CHK_1_H["sample"] = "CHK_1_H"
IE_1_H["sample"] = "IE_1_H"

# CHK_2_H = pd.read_csv("data/chicken_galGal6/CHK_2_H_merged_with_gene.csv")
IE_2_H = pd.read_csv("data/IE_galGal6/IE_2_H_merged_with_gene.csv")

# CHK_2_H["sample"] = "CHK_2_H"
IE_2_H["sample"] = "IE_2_H"

# CHK_3_H = pd.read_csv("data/chicken_galGal6/CHK_3_H_merged_with_gene.csv")
# CHK_4_H = pd.read_csv("data/chicken_galGal6/CHK_4_H_merged_with_gene.csv")
# CHK_5_H = pd.read_csv("data/chicken_galGal6/CHK_5_H_merged_with_gene.csv")

# CHK_3_H["sample"] = "CHK_3_H"
# CHK_4_H["sample"] = "CHK_4_H"
# CHK_5_H["sample"] = "CHK_5_H"

IE_1_H["Animal"] = "eagle"
IE_2_H["Animal"] = "eagle"

CHK_1_H["Animal"] = "chicken"
# CHK_2_H["Animal"] = "chicken"
# CHK_3_H["Animal"] = "chicken"
# CHK_4_H["Animal"] = "chicken"
# CHK_5_H["Animal"] = "chicken"

frames = [IE_1_H, CHK_1_H, IE_2_H]

merged = pd.concat(frames, ignore_index=True)

merged.to_csv("data/chick_comp/chick_eagle_merged_with_gene_eag.csv", index=False)

animal1 = "chicken"
animal2 = "eagle"

sites = merged[merged["methyl_status"].str.isupper()]

cols = ["chromosome", "position", "strand", "methyl_status"]

a = set(
    sites.loc[sites["Animal"] == animal1, cols]
      .itertuples(index=False, name=None)
)

b = set(
    sites.loc[sites["Animal"] == animal2, cols]
      .itertuples(index=False, name=None)
)

shared = len(a & b)
only_a = len(a - b)
only_b = len(b - a)

print(f"Shared: {shared}")
print(f"Only {animal1}: {only_a}")
print(f"Only {animal2}: {only_b}")

union = len(a | b)

jaccard = shared / union

print("Jaccard similarity =", jaccard)

fasta_file = "D:/bird_methylome/CHK/galGal6_genome.fasta"

total_c = 0

for record in SeqIO.parse(fasta_file, "fasta"):
    total_c += str(record.seq).upper().count("C")

print("Total cytosines:", total_c)

A = len(a)
B = len(b)

pval = hypergeom.sf(shared - 1, total_c, A, B)

print("Hypergeometric p-value:", f"{pval:.4f}")

# expected overlap
expected = (A * B) / total_c

# fold enrichment
fold_enrichment = shared / expected

print(f"Observed overlap: {shared:,}")
print(f"Expected overlap: {expected:.2f}")
print(f"Fold enrichment: {fold_enrichment:.2f}")

# now to do within eagle comparisons
cols = ["chromosome", "position", "strand", "methyl_status"]

IE_1 = "IE_1_H"
IE_2 = "IE_2_H"

a_IE = set(
    sites.loc[sites["sample"] == IE_1, cols]
      .itertuples(index=False, name=None)
)

b_IE = set(
    sites.loc[sites["sample"] == IE_2, cols]
      .itertuples(index=False, name=None)
)

shared_IE = len(a_IE & b_IE)
only_a_IE = len(a_IE - b_IE)
only_b_IE = len(b_IE - a_IE)

print(f"Shared: {shared_IE}")
print(f"Only {IE_1}: {only_a_IE}")
print(f"Only {IE_2}: {only_b_IE}")

union_IE = len(a_IE | b_IE)

jaccard_IE = shared_IE / union_IE

print("Jaccard similarity =", jaccard_IE)

A_IE = len(a_IE)
B_IE = len(b_IE)

pval_IE = hypergeom.sf(shared_IE - 1, total_c, A_IE, B_IE)

print("Hypergeometric p-value:", f"{pval_IE:.4f}")

# expected overlap
expected_IE = (A_IE * B_IE) / total_c

# fold enrichment
fold_enrichment_IE = shared_IE / expected_IE

print(f"Observed overlap: {shared_IE:,}")
print(f"Expected overlap: {expected_IE:.2f}")
print(f"Fold enrichment: {fold_enrichment_IE:.2f}")
