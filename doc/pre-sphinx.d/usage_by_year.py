#!/usr/bin/env python
"""Create barplot of AMICI publications by year."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main():
    script_path = Path(__file__).parent
    outfile = script_path.parent / "gfx" / "usage_by_year.png"

    # set rcParams for better readability
    plt.rcParams.update(
        {
            "figure.figsize": (5, 3),
            "figure.dpi": 150,
            "axes.titlesize": 14,
            "axes.labelsize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
        }
    )

    df = pd.read_csv(script_path.parent / "usage_by_year.csv")
    plt.bar(df.year, df.citations)
    plt.xlabel("Year")
    plt.ylabel("Citations")
    plt.tight_layout()
    plt.savefig(outfile)


if __name__ == "__main__":
    main()
