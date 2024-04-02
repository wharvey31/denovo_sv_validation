#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: William T. Harvey

import pandas as pd


if __name__ == "__main__":

    bed_df = pd.read_csv(snakemake.input.bed, sep="\t")
    for i, file in enumerate(snakemake.input.val):
        if i == 0:
            df = pd.read_csv(file, sep="\t")
        else:
            df = df.merge(pd.read_csv(file, sep="\t"), how="outer")

    df["DNSVAL"] = df.apply(
        lambda row: "PASS"
        if list(set([row[x] for x in df.columns if x != "ID"])) == ["VALID"]
        else "FAIL",
        axis=1,
    )

    df.loc[df['DNSVAL'] == "PASS"][["ID", "DNSVAL"]].merge(bed_df).to_csv(
        snakemake.output.val, sep="\t", index=False
    )
    df.to_csv(snakemake.output.all_calls, sep="\t", index=False)
