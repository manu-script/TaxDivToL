import pandas as pd
import re
from joblib import Parallel, delayed
import os
from functools import reduce
import sys
from multiprocessing import cpu_count

mmseqs = sys.argv[1]


def lca_process(lca):
    part = lca.split('/')[-1].split('_')[0]
    df = pd.read_csv(lca, sep='\t', header=None, names=['query', 'lineage'], usecols=[0, 4])
    df.lineage = df.lineage.map(lambda x: ';'.join(re.findall(r'[dkpcofgs]_[A-Z][a-z]+[^;]*', x)))
    df = df[df.lineage != '']
    dup = df[df.duplicated(subset='query', keep=False)].copy()
    df = df.loc[~df.duplicated(subset='query', keep=False)]
    df.drop('query', axis=1, inplace=True)
    dup.lineage = dup.lineage.map(lambda x: [''.join(re.findall(r'(d_[A-Z][a-z]+[^;]+)', x))] + [''.join(re.findall(r'(k_[A-Z][a-z]+[^;]+)', x))] +
                                            [''.join(re.findall(r'(p_[A-Z][a-z]+[^;]+)', x))] + [''.join(re.findall(r'(c_[A-Z][a-z]+[^;]+)', x))] +
                                            [''.join(re.findall(r'(o_[A-Z][a-z]+[^;]+)', x))] + [''.join(re.findall(r'(f_[A-Z][a-z]+[^;]+)', x))] +
                                            [''.join(re.findall(r'(g_[A-Z][a-z]+[^;]+)', x))] + [''.join(re.findall(r'(s_[A-Z][a-z]+[^;]+)', x))])
    dup = pd.DataFrame(map(lambda x: ';'.join([a for a, b in zip(x[1].iloc[0, 1], x[1].iloc[1, 1]) if a == b != '']), dup.groupby('query')), columns=['lineage'])
    df = df.append(dup, ignore_index=True)
    df = df[df.lineage != '']
    df = df.lineage.map(lambda x: x.split(';')).map(lambda y: [';'.join(y[:len(y) - i]) for i in range(len(y))]).apply(pd.Series).stack().reset_index(
        drop=True).value_counts().rename_axis('lineage').reset_index(name=part)
    return df[['lineage', part]]


parts = Parallel(n_jobs=cpu_count())(delayed(lca_process)(l) for l in [os.path.join(mmseqs, i) for i in os.listdir(mmseqs) if i.endswith('lca.tsv.gz')])
parts = reduce(lambda left, right: pd.merge(left, right, on=['lineage'], how='outer'), parts).fillna(0)
parts.loc[:, parts.columns != 'lineage'] = parts.loc[:, parts.columns != 'lineage'].astype(int)
parts.to_csv(os.path.join(mmseqs, 'counts.tsv'), sep='\t', index=False)
