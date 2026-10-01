import pandas as pd
import numpy as np

CpG = pd.read_csv("D:/bird_methylome/IE/IE_galGal6/IE_1_H_CpG_context_file.txt", sep='\t')
CHH = pd.read_csv("D:/bird_methylome/IE/IE_galGal6/IE_1_H_CHH_context_file.txt", sep='\t')
CHG = pd.read_csv("D:/bird_methylome/IE/IE_galGal6/IE_1_H_CHG_context_file.txt", sep='\t')
cov = pd.read_csv("D:/bird_methylome/IE/IE_galGal6/IE_1_H_coverage_for_all.tabular", sep='\t')

CpG.columns = ['read_index', 'strand', 'chromosome', 'position', 'methyl_status']
CHH.columns = ['read_index', 'strand', 'chromosome', 'position', 'methyl_status']
CHG.columns = ['read_index', 'strand', 'chromosome', 'position', 'methyl_status']
cov.columns = ['chromosome', 'position', 'position_end', 'coverage_percent', 'methylated_reads', 'unmethylated_reads']

print("number of methylated sites in unfiltered coverage: ", len(cov))

cov2 = cov[(cov['methylated_reads'] + cov['unmethylated_reads'] >= 10) & (cov['methylated_reads'] > 10)]

print("number of methylated sites in filtered coverage: ", len(cov2))

df_merged1 = pd.merge(cov2, CpG, CHH, CHG, on=["chromosome", "position", "strand"], how="outer")

df_merged1.to_csv("data/IE_galGal6/IE_1_H_merged.csv")