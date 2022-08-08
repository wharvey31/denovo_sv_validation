
rule int_check:
    input:
        bed=find_bed,
        int_val=find_int,
    output:
        tab="temp/validation/SVPOP/{val_type}/{vartype}_{svtype}/{sample}_raw.tsv",
    resources:
        mem=4,
        hrs=24,
    threads: 1
    run:
        val_dict = {"A": "VALID", "A,B": "NOTVALID"}
        val_df = pd.merge(
            pd.read_csv(
                input.int_val[0], sep="\t", usecols=["ID_A", "SOURCE_SET"]
            ).dropna(),
            pd.read_csv(
                input.int_val[1], sep="\t", usecols=["ID_A", "SOURCE_SET"]
            ).dropna(),
            on="ID_A",
            suffixes=["_MO", "_FA"],
        )
        val_df = val_df.rename(columns={"ID_A": "ID"})
        all_df = pd.read_csv(input.bed, sep="\t").merge(val_df)
        val_df[f"SVPOP_{wildcards.val_type}_MO"] = val_df.apply(
            lambda row: val_dict[row["SOURCE_SET_MO"]], axis=1
        )
        val_df[f"SVPOP_{wildcards.val_type}_FA"] = val_df.apply(
            lambda row: val_dict[row["SOURCE_SET_FA"]], axis=1
        )
        val_df[
            ["ID", f"SVPOP_{wildcards.val_type}_MO", f"SVPOP_{wildcards.val_type}_FA"]
        ].to_csv(output.tab, sep="\t", index=False)
