#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: William T. Harvey


if __name__ == "__main__":

    for i, file in enumerate(snakemake.input.val):
        if i == 0:
            df = pd.read_csv(file, sep="\t")
        else:
            df = df.merge(pd.read_csv(file, sep="\t"))
    sample_cols = [
        x for x in df.columns if snakemake.wildcards.sample in x and x != "ID"
    ]
    parental_cols = [
        x for x in df.columns if snakemake.wildcards.sample not in x and x != "ID"
    ]
    sample_val = df[sample_cols].str.replace("PASS", True)
    parent_val = df[parental_cols].str.replace("PASS", True)
    df.to_csv(snakemake.output.all_calls, sep="\t", index=False)
    sample_val["ID"].merge(parent_val["ID"]).to_csv(
        snakemake.output.val, sep="\t", index=False
    )
