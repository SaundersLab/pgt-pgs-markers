import glob
import json
import os
import pandas as pd

rows = []
for path in glob.glob('reports/alignment/*.json'):
    with open(path) as f:
        sample_dict = json.load(f)['QC-passed reads']
        sample_dict['Sample Accession'] = os.path.splitext(
            os.path.basename(path))[0]
        rows.append(sample_dict)
df = pd.DataFrame(rows)

metadata = pd.read_csv('data/metadata/metadata.csv')
metadata['fsp'] = metadata.Host.apply(
    lambda s: 'Pgt' if 'heat' in s else {'Oat': 'Pga', 'Rye': 'Pgs'}[s]
)
df = pd.merge(metadata, df, on='Sample Accession')

groupby = df.groupby('fsp')
pd.DataFrame({
    'mean': groupby['mapped %'].mean(),
    'std': groupby['mapped %'].std(),
}).to_csv('reports/mapping_summary.csv', index=None)
df = df[
    [
        'Sample Accession',
        'Sample Name',
        'paired in sequencing',
        'mapped',
        'mapped %'
    ]
]
df = df.rename(
    columns={
        'paired in sequencing': 'Total number of reads',
        'mapped': 'Number of mapped reads',
        'mapped %': 'Mapped reads (%)'
    }
)
df.to_csv('reports/mapping.csv', index=None)
