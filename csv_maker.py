import pandas as pd
import numpy as np

CpG = pd.read_csv("D:/bird_methylome/CHK/chick_galGal6/CpG_context_chicken.txt", sep='\t')
CHH = pd.read_csv("D:/bird_methylome/CHK/chick_galGal6/CHH_context_chicken.txt", sep='\t')
CHG = pd.read_csv("D:/bird_methylome/CHK/chick_galGal6/CHG_context_chicken.txt", sep='\t')
cov = pd.read_csv("D:/bird_methylome/CHK/chick_galGal6/chicken_coverage_for_all.tabular", sep='\t')

CpG.columns = ['read_index', 'strand', 'chromosome', 'position', 'methyl_status']
CHH.columns = ['read_index', 'strand', 'chromosome', 'position', 'methyl_status']
CHG.columns = ['read_index', 'strand', 'chromosome', 'position', 'methyl_status']
cov.columns = ['chromosome', 'position', 'position_end', 'coverage_percent', 'methylated_reads', 'unmethylated_reads']

print("number of methylated sites in unfiltered coverage: ", len(cov))

cov2 = cov[(cov['methylated_reads'] + cov['unmethylated_reads'] >= 10) & (cov['methylated_reads'] > 10)]

print("number of methylated sites in filtered coverage: ", len(cov2))

df_merged1 = pd.merge(CpG, cov2, on=["chromosome", "position"], how="outer")

df_merged1 = df_merged1.drop(columns=["read_index"])

df_merged2 = pd.merge(CHH, df_merged1, on=["chromosome", "position", "strand", "methyl_status"], how="outer")

df_merged2 = df_merged2.drop(columns=["read_index"])

df_merged3 = pd.merge(CHG, df_merged2, on=["chromosome", "position", "strand", "methyl_status"], how="outer")

df_merged4 = df_merged3.dropna()

df_merged5 = df_merged4.drop(columns=["read_index"])

df_merged6 = df_merged5.drop_duplicates(subset=["chromosome", "position", "strand", "methyl_status"], keep="first")

df_merged6.to_csv("data/chicken_galGal6/chicken_galGal6_merged.csv", index=False)