import pandas as pd
import sys

coverage_table_path = sys.argv[1]
coverage_stats_table_path = sys.argv[2]

pgs_samples = ['SAMEA9973400', 'SAMEA9973401', 'SAMEA9973402', 'SAMEA9973403']
pga_samples = ['SAMEA104219957', 'SAMEA104219958']

coverage_table = pd.read_csv(coverage_table_path).drop(columns=pga_samples, errors='ignore')

pgt_samples = sorted(set(coverage_table) - set(pgs_samples))
pd.DataFrame({
    'pgt_min': coverage_table[pgt_samples].min(axis=1),
    'pgs_min': coverage_table[pgs_samples].min(axis=1),
    'pgs_max': coverage_table[pgs_samples].max(axis=1),
}).to_csv(coverage_stats_table_path, index=None)
