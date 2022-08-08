#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: William T. Harvey

import pandas as pd
from pybedtools import BedTool


if __name__ == "__main__":

    df = pd.read_csv(snakemake.input.tab, sep="\t")
    call = BedTool(snakemake.input.bed)
    regions = (
        BedTool.from_dataframe(df[["#CHROM", "POS", "END", "ID"]])
        .intersect(call, wa=True)
        .to_dataframe(names=["#CHROM", "POS", "END", "ID"])
    )
    regions[
        f"CALLABLE_{snakemake.wildcards.val_type}_{snakemake.wildcards.parent}_{snakemake.wildcards.hap}"
    ] = "VALID"
    df = df.merge(regions, how="left")
    df[
        f"CALLABLE_{snakemake.wildcards.val_type}_{snakemake.wildcards.parent}_{snakemake.wildcards.hap}"
    ] = df[
        f"CALLABLE_{snakemake.wildcards.val_type}_{snakemake.wildcards.parent}_{snakemake.wildcards.hap}"
    ].fillna(
        "NOTVALID"
    )

    df[
        [
            "ID",
            f"CALLABLE_{snakemake.wildcards.val_type}_{snakemake.wildcards.parent}_{snakemake.wildcards.hap}",
        ]
    ].to_csv(snakemake.output.tab, sep="\t", index=False)
