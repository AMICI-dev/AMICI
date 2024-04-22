#!/usr/bin/env python3

"""
Aggregate computation times from different benchmarks and plot
"""

import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# read benchmark results for different models

outfile = "computation_times.csv"
df = pd.concat(
    [
        pd.read_csv(f, header=[0], index_col=[0])
        .rename(columns={"0": "_".join(f.split("_")[:2])})
        .T
        for f in os.listdir()
        if f.endswith(".csv")
        if f != outfile
    ]
)
df.sort_values("np", inplace=True)

df.to_csv(outfile)


ratios = (
    pd.concat(
        [df[sensi] / df["t_sim"].values for sensi in ["t_fwd", "t_adj"]]
        + [df.np],
        axis=1,
    )
    .reset_index()
    .melt(id_vars=["index", "np"])
    .rename(
        columns={"index": "model", "variable": "sensitivity", "value": "ratio"}
    )
)
ratios["sensitivity"] = ratios["sensitivity"].replace(
    {"t_fwd": "forward", "t_adj": "adjoint"}
)


plt.figure(figsize=(10, 5))
g = sns.barplot(
    data=ratios, order=list(df.index), x="model", y="ratio", hue="sensitivity"
)
for ir, row in ratios.iterrows():
    if row.sensitivity == "adjoint":
        continue
    g.text(
        ir,
        row["np"],
        int(row["np"]),
        color="black",
        ha="center",
        weight="bold",
    )

plt.xticks(rotation=30, horizontalalignment="right")
plt.tight_layout()
plt.savefig("computation_times.png")
