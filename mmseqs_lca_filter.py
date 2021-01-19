import pandas as pd
import re
import numpy as np
import os
from functools import reduce
import sys
import subprocess


def make_krona(df):
    df.lineage = df.lineage.map(lambda x: [''.join(re.findall(r'(d_[A-Z][a-z]+[^;]+)', x))] + [''.join(re.findall(r'(k_[A-Z][a-z]+[^;]+)', x))] +
                                          [''.join(re.findall(r'(p_[A-Z][a-z]+[^;]+)', x))] + [''.join(re.findall(r'(c_[A-Z][a-z]+[^;]+)', x))] +
                                          [''.join(re.findall(r'(o_[A-Z][a-z]+[^;]+)', x))] + [''.join(re.findall(r'(f_[A-Z][a-z]+[^;]+)', x))] +
                                          [''.join(re.findall(r'(g_[A-Z][a-z]+[^;]+)', x))] + [''.join(re.findall(r'(s_[A-Z][a-z]+[^;]+)', x))])
    df[['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']] = pd.DataFrame(df.lineage.tolist(), index=df.index)
    df = df[df['genus'] != '']
    df[['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus']].to_csv('krona.txt', sep='\t', header=False, index=False)
    subprocess.run(['ktImportText', '-n', 'Tree of Life', '-o', 'krona.html', '-q', 'krona.txt'], check=True)


mmseqs = [pd.read_csv(os.path.join(i, 'counts.tsv'), sep='\t', header=0) for i in sys.argv[1:]]
union = reduce(lambda left, right: pd.merge(left, right, on=['lineage'], how='outer'), mmseqs).fillna(0)
union['count'] = union.loc[:, union.columns != 'lineage'].sum(axis=1)
union = union[['lineage', 'count']]
union['parent'] = union.lineage.map(lambda x: union['count'][union.lineage == ';'.join(x.split(';')[:-1])]).map(lambda y: y.values[0] if len(y) == 1 else np.nan)
filtered = union.lineage[(union.lineage.map(lambda x: bool(re.search(r';f_[A-Z][a-z]+', x)) or bool(re.search(r';g_[A-Z][a-z]+', x)) or bool(re.search(r';s_[A-Z][a-z]+', x)))) & ((union['count'] <= 10) | (union['count'] < union['parent'] * 0.01))]
union.drop(filtered.map(lambda x: union.index[(union.lineage == x) | (union.lineage.str.startswith(x + ';'))]).apply(pd.Series).stack().reset_index(drop=True).astype(int).unique(), inplace=True)
union = union[['lineage', 'count']]
union['count'] = union['count'].astype(int)
union.to_csv('counts.tsv', sep='\t', header=True, index=False)
make_krona(union.copy())
union[union.lineage.map(lambda x: bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Tree_of_life_genus_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Bacteria;') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Bacteria_genus_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Viruses;') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Viruses_genus_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Archaea;') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Archaea_genus_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Eukaryota_genus_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;') and not x.startswith('d_Eukaryota;k_') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Protists_genus_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Metazoa_genus_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Fungi;') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Fungi_genus_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Viridiplantae;') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Viridiplantae_genus_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Chordata_genus_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;') and not x.startswith('d_Eukaryota;k_Metazoa;p_Chordata') and bool(
    re.search(r';g_[A-Z][a-z]+', x)) and not bool(re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Invertebrates_genus_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Arthropoda') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Arthropoda_genus_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;c_Mammalia') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Mammalia_genus_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;c_Amphibia') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Amphibia_genus_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;c_Aves') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Aves_genus_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;c_Actinopteri') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Actinopteri_genus_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;c_Chondrichthyes') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Chondrichthyes_genus_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: (x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;o_Crocodylia') or x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;o_Lepidosauria') or
                                   x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;o_Testudines')) and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Reptiles_genus_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Tree_of_life_family_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Bacteria;') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Bacteria_family_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Viruses;') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Viruses_family_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Archaea;') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Archaea_family_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Eukaryota_family_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;') and not x.startswith('d_Eukaryota;k_') and bool(
    re.search(r';f_[A-Z][a-z]+', x)) and not bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Protists_family_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(
    re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Metazoa_family_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Fungi;') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(
    re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Fungi_family_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Viridiplantae;') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(
    re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Viridiplantae_family_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(
    re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Chordata_family_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;') and not x.startswith('d_Eukaryota;k_Metazoa;p_Chordata') and bool(
    re.search(r';g_[A-Z][a-z]+', x)) and not bool(re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Invertebrates_family_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Arthropoda') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(
    re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Arthropoda_family_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;c_Mammalia') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(
    re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Mammalia_family_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;c_Amphibia') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(
    re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Amphibia_family_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;c_Aves') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(
    re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Aves_family_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;c_Actinopteri') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(
    re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Actinopteri_family_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;c_Chondrichthyes') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(
    re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Chondrichthyes_family_level.tsv', sep='\t', header=True, index=False)

union[union.lineage.map(lambda x: (x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;o_Crocodylia') or x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;o_Lepidosauria') or
                                   x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;o_Testudines')) and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(
    re.search(r';g_[A-Z][a-z]+', x)) and not bool(
    re.search(r';s_[A-Z][a-z]+', x)))].to_csv('Reptiles_family_level.tsv', sep='\t', header=True, index=False)


for i, parts in enumerate(mmseqs):
    parts = parts[parts.lineage.isin(union.lineage)]
    parts.loc[:, parts.columns != 'lineage'] = parts.loc[:, parts.columns != 'lineage'].applymap(lambda x: 1 if x != 0 else 0)
    parts[parts.lineage.map(lambda x: bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Tree_of_life_genus_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Bacteria;') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Bacteria_genus_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Viruses;') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Viruses_genus_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Archaea;') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Archaea_genus_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Eukaryota_genus_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;') and not x.startswith('d_Eukaryota;k_') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Protists_genus_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Metazoa_genus_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Fungi;') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Fungi_genus_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Viridiplantae;') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Viridiplantae_genus_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Chordata_genus_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;') and not x.startswith('d_Eukaryota;k_Metazoa;p_Chordata') and bool(
        re.search(r';g_[A-Z][a-z]+', x)) and not bool(re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Invertebrates_genus_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Arthropoda') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Arthropoda_genus_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;c_Mammalia') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Mammalia_genus_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;c_Amphibia') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Amphibia_genus_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;c_Aves') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Aves_genus_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;c_Actinopteri') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Actinopteri_genus_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;c_Chondrichthyes') and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Chondrichthyes_genus_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: (x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;o_Crocodylia') or x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;o_Lepidosauria') or
                                       x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;o_Testudines')) and bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Reptiles_genus_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Tree_of_life_family_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Bacteria;') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Bacteria_family_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Viruses;') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Viruses_family_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Archaea;') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Archaea_family_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Eukaryota_family_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;') and not x.startswith('d_Eukaryota;k_') and bool(
        re.search(r';f_[A-Z][a-z]+', x)) and not bool(re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Protists_family_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(
        re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Metazoa_family_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Fungi;') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(
        re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Fungi_family_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Viridiplantae;') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(
        re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Viridiplantae_family_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(
        re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Chordata_family_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;') and not x.startswith('d_Eukaryota;k_Metazoa;p_Chordata') and bool(
        re.search(r';g_[A-Z][a-z]+', x)) and not bool(re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Invertebrates_family_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Arthropoda') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(
        re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Arthropoda_family_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;c_Mammalia') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(
        re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Mammalia_family_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;c_Amphibia') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(
        re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Amphibia_family_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;c_Aves') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(
        re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Aves_family_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;c_Actinopteri') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(
        re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Actinopteri_family_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;c_Chondrichthyes') and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(
        re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Chondrichthyes_family_level.tsv'), sep='\t', header=True, index=False)

    parts[parts.lineage.map(lambda x: (x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;o_Crocodylia') or x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;o_Lepidosauria') or
                                       x.startswith('d_Eukaryota;k_Metazoa;p_Chordata;o_Testudines')) and bool(re.search(r';f_[A-Z][a-z]+', x)) and not bool(
        re.search(r';g_[A-Z][a-z]+', x)) and not bool(
        re.search(r';s_[A-Z][a-z]+', x)))].to_csv(os.path.join(sys.argv[i+1], 'Reptiles_family_level.tsv'), sep='\t', header=True, index=False)

