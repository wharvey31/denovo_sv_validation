# subseq_table_sv_bed
#
# Make SV validation table.
#
# Wildcards:
# * sample: Biological sample
# * set_def: Entry in SET_DEF (without the sv prefix)
# * caller: Variant caller
# * svtype: ins or del
rule subseq_table_setdef_bed:
    input:
        bed=find_bed,
    output:
        bed="temp/tables/subseq/{sample}/{parent}_{val_type}/{set_def}/setdef_{vartype}_{svtype}.bed.gz",
    resources:
        mem=4,
        hrs=24,
    threads: 1
    run:
        # Get parameters
        if wildcards.set_def not in SET_DEF:
            raise RuntimeError(
                "Set definition not found in SET_DEF: " + wildcards.set_def
            )

        svlen_min, svlen_max, win_size = SET_DEF[wildcards.set_def]

        # Process variants
        if "bed" in input.keys() and bool(
            input.bed
        ):  # Function returns an empty list if no input, and Snakemake omits those from the input object

            df = pd.read_csv(
                input.bed,
                sep="\t",
                usecols=("#CHROM", "POS", "END", "ID", "SVTYPE", "SVLEN"),
            )

            df_fai = get_ref_fai(REF_FA + ".fai")

            # Filter chromosomes that were not aligned to (alt and decoys for short-read callers)
            df = df.loc[df["#CHROM"].apply(lambda val: val in df_fai.index)]

            # Subset
            if svlen_min is not None:
                df = df.loc[df["SVLEN"] >= svlen_min]

            if svlen_max is not None:
                df = df.loc[df["SVLEN"] < svlen_max]

            # Add sample and caller
            df["SAMPLE"] = wildcards.sample
            df["CALLER"] = "POTENTIAL_DENOVO"

            # Annotate window (replace POS and END with the window coordinates)
            df["WIN_FLANK"] = win_size

            # Define windows, replace POS and END
            df["ORG_POS"] = df["POS"]
            df["ORG_END"] = df["END"]

            if df.shape[0] > 0:
                df["POS"] = df.apply(
                    lambda row: np.max([row["POS"] - row["WIN_FLANK"], 0]), axis=1
                )
                df["END"] = df.apply(
                    lambda row: np.min(
                        [row["END"] + row["WIN_FLANK"], df_fai[row["#CHROM"]]]
                    ),
                    axis=1,
                )

            df["WIN_L"] = df["ORG_POS"] - df["POS"]
            df["WIN_R"] = df["END"] - df["ORG_END"]

            df["WINDOW_SIZE"] = df["WIN_L"] + df["WIN_R"]

            del (df["ORG_POS"], df["ORG_END"])

            if df.shape[0] > 0:
                df["WINDOW"] = df.apply(
                    lambda row: "{#CHROM}:{POS}-{END}".format(**row), axis=1
                )
            else:
                df["WINDOW"] = []
        else:
            # No variant calls

            df = pd.DataFrame(
                [],
                columns=[
                    "#CHROM",
                    "POS",
                    "END",
                    "ID",
                    "SVTYPE",
                    "SVLEN",
                    "SAMPLE",
                    "CALLER",
                    "WIN_FLANK",
                    "WIN_L",
                    "WIN_R",
                    "WINDOW_SIZE",
                    "WINDOW",
                ],
            )

        # Write
        df.to_csv(output.bed, sep="\t", index=False, compression="gzip")

# subseq_tab_window_single
#
# Get stats for a single window.
rule subseq_tab_window_single:
    input:
        bed=rules.subseq_table_setdef_bed.output.bed,
        aln=lambda wildcards: get_aln_source(wildcards, config["READS"]),
    output:
        tsv=temp(
            "temp/tables/sample/{sample}/{parent}_{val_type}/{set_def}/{vartype}_{svtype}/{parent}_{val_type}.tsv.gz"
        ),
    resources:
        mem=4,
        hrs=24,
    threads: 1
    singularity:
        "docker://eichlerlab/subseqpy:1.0"
    script:
        "../scripts/subseq_subp.py"


#
# Subsequence tables
#


# subseq_merge_setdef
#
# Merge subseq results across set definitions (one table for variant callset vs one alignment source).
rule subseq_merge_setdef:
    input:
        tsv=gather_setdef,
    output:
        tsv="temp/tables/subseq/{sample}/{parent}_{val_type}/{vartype}_{svtype}/{parent}_{val_type}.tsv.gz",
    resources:
        mem=4,
        hrs=24,
    threads: 1
    run:
        return pd.concat(
            [pd.read_csv(tsv_file_name, sep="\t") for tsv_file_name in input.tsv],
            axis=0,
        ).to_csv(output.tsv, sep="\t", index=False, compression="gzip")


#
# Validation tables
#


# subseq_val
#
# Make a validation table
rule subseq_val:
    input:
        tsv=rules.subseq_merge_setdef.output.tsv,
    output:
        tsv="temp/tables/validation/{sample}/{parent}_{val_type}/{vartype}_{svtype}.tsv.gz",
    resources:
        mem=4,
        hrs=24,
    threads: 1
    run:
        # # Check alignsource generator function
        # if ALNSOURCE_PLOIDY_DICT[wildcards.alnsource] != align_summary_diploid:
        #     raise RuntimeError('Validation with size50 requires results from align_summary_diploid')

        # Read
        df = pd.read_csv(input.tsv, sep="\t")

        if wildcards.vartype == "sv":
            strategy = "size50_2_4"
        else:
            strategy = "size20_2_4"

        # Validate
        df = validate_summary(df, strategy=strategy)

        # Write
        df.to_csv(output.tsv, sep="\t", index=False, compression="gzip")


rule gather_parents:
    input:
        fa_val=subseq_father,
        mo_val=subseq_mother,
    output:
        subseq="temp/validation/READS/{val_type}/{vartype}_{svtype}/{sample}_raw.tsv",
    resources:
        mem=4,
        hrs=24,
    threads: 1
    run:
        df = pd.merge(
            pd.read_csv(input.fa_val, sep="\t"),
            pd.read_csv(input.mo_val, sep="\t"),
            on="ID",
            suffixes=[f"_{wildcards.val_type}_fa", f"_{wildcards.val_type}_mo"],
        )
        # Swap validation structure
        val_columns = []
        for parent in ["mo", "fa"]:
            df[f"VAL_{wildcards.val_type}_{parent}"] = df[
                f"VAL_{wildcards.val_type}_{parent}"
            ].replace({"VALID": "NOTVALID", "NOTVALID": "VALID"})
            val_columns.append(f"VAL_{wildcards.val_type}_{parent}")
        df[["ID"] + val_columns].to_csv(output.subseq, sep="\t", index=False)
