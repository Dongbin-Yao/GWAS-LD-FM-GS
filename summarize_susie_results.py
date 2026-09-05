#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
summarize_susie_results.py
==========================

Summarize and rank SNPs according to SuSiE posterior inclusion
probabilities (PIP).

Input
-----
One or more SuSiE PIP result files containing:

    SNP
    PIP

Output
------
A combined SNP-level PIP ranking.

For duplicated SNPs, the maximum PIP is retained.
"""

import argparse
from pathlib import Path

import pandas as pd


def read_pip_file(filename):

    filename = Path(filename)

    if not filename.exists():
        raise FileNotFoundError(
            f"File not found: {filename}"
        )

    df = pd.read_csv(
        filename,
        sep=r"\s+",
        engine="python"
    )

    required = {"SNP", "PIP"}

    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            f"{filename} is missing columns: "
            f"{', '.join(sorted(missing))}"
        )

    df["PIP"] = pd.to_numeric(
        df["PIP"],
        errors="coerce"
    )

    df = df.dropna(
        subset=["SNP", "PIP"]
    )

    return df[["SNP", "PIP"]]


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Combine and rank SuSiE PIP results."
        )
    )

    parser.add_argument(
        "--input",
        nargs="+",
        required=True,
        help="One or more SuSiE PIP result files."
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output PIP ranking file."
    )

    args = parser.parse_args()

    all_results = []

    for filename in args.input:

        df = read_pip_file(filename)

        all_results.append(df)

    combined = pd.concat(
        all_results,
        ignore_index=True
    )

    # If a SNP appears more than once, retain the highest PIP.
    combined = (
        combined
        .groupby("SNP", as_index=False)["PIP"]
        .max()
    )

    combined = combined.sort_values(
        ["PIP", "SNP"],
        ascending=[False, True]
    )

    combined["PIP_Rank"] = range(
        1,
        len(combined) + 1
    )

    combined = combined[
        ["PIP_Rank", "SNP", "PIP"]
    ]

    output = Path(args.output)

    output.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    combined.to_csv(
        output,
        sep="\t",
        index=False
    )

    print(
        f"Combined SNPs: {len(combined)}"
    )

    print(
        f"Output: {output}"
    )


if __name__ == "__main__":
    main()
