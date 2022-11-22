import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# read benchmark results for different models
df = pd.concat([
    pd.read_csv(f, header=[0], index_col=[0]).rename(columns={'0': '_'.join(f.split('_')[:2])}).T
    for f in os.listdir() if f.endswith('.csv')
])

# regex applied to README.md to extract model sizes
# \| \[([\w_]+)\]\([\w\-\/]+\)\s+\|\s+[0-9]+\s\|\s+[0-9]+\s\|\s+0\s\|\s+[0-9]\s\|\s+[0-9]\s\|\s+[0-9]+\s\|\s+[0-9]+\s\| normal\s+\|\s+([0-9]+) \| \[\\\[1\\\]\]\(http://identifiers.org/(doi|pubmed)/[0-9\./\w-]+\)\s+\|
# '\1': \2,
n_species = {
    'Alkan_SciSignal2018': 36,
    'Bachmann_MSB2011': 25,
    'Beer_MolBioSystems2014': 4,
    'Bertozzi_PNAS2020': 3,
    'Blasi_CellSystems2016': 16,
    'Boehm_JProteomeRes2014': 8,
    'Borghans_BiophysChem1997': 3,
    'Brannmark_JBC2010': 9,
    'Bruno_JExpBot2016': 7,
    'Chen_MSB2009': 500,
    'Crauste_CellSystems2017': 5,
    'Elowitz_Nature2000': 8,
    'Fiedler_BMC2016': 6,
    'Froehlich_CellSystems2018': 1396,
    'Fujita_SciSignal2010': 9,
    'Giordano_Nature2020': 13,
    'Isensee_JCB2018': 25,
    'Laske_PLOSComputBiol2019': 41,
    'Lucarelli_CellSystems2018': 33,
    'Okuonghae_ChaosSolitonsFractals2020': 9,
    'Oliveira_NatCommun2021': 9,
    'Perelson_Science1996': 4,
    'Rahman_MBS2016': 7,
    'Raimundez_PCB2020': 22,
    'SalazarCavazos_MBoC2020': 75,
    'Schwen_PONE2014': 11,
    'Sneyd_PNAS2002': 6,
    'Weber_BMC2015': 7,
    'Zhao_QuantBiol2020': 5,
    'Zheng_PNAS2012': 15,
}

df = df.apply(np.log10)
df = df.sub(df.AMICI.values, axis=0)

plt.figure(figsize=(10, 5))
g = sns.heatmap(
    # sort according to # species
    df.loc[sorted(df.index, key=lambda x: n_species[x]), :],
    fmt='.1f', cmap='RdBu_r', vmin=-2, vmax=2, annot=True,
    cbar_kws={'label': 'log10 foldchange simulation time'},
)
plt.tight_layout()
plt.show()
plt.savefig('benchmark-models-julia.pdf')
