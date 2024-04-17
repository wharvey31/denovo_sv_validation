#!/bin/env python3

import os
import pandas as pd
import io


def get_len_list(window, aln_input, subseq_exe):
    """
    Get a list of alignment record lengths over a window.

    :param window: Position string (chrom:pos-end). Coordinates are 1-based inclusive (not BED).
    :param aln_input: Alignment input as BAM or CRAM.
    :param subseq_exe: Path to subseq executable.
    """

    # Open process
    proc = subprocess.Popen(
        [subseq_exe, "-b", "-r", window, aln_input], stdout=subprocess.PIPE
    )

    stdout, stderr = proc.communicate()

    # Return list
    with io.StringIO(stdout.decode()) as in_file:
        return [len(record.seq) for record in SeqIO.parse(in_file, "fasta")]


# Read BED



# Read region BED
region_bed = pd.read_csv(snakemake.input.bed, sep="\t")

region_bed["ALNSAMPLE"] = snakemake.wildcards.parent
region_bed["ALNSOURCE"] = snakemake.wildcards.val_type

region_bed = region_bed.loc[
    :,
    [
        "ID",
        "SVTYPE",
        "SVLEN",
        "WINDOW",
        "WINDOW_SIZE",
        "SAMPLE",
        "CALLER",
        "ALNSAMPLE",
        "ALNSOURCE",
    ],
].reset_index(drop=True)

# Get alignment source
aln_input_exists = os.path.isfile(snakemake.input.aln)

# Get summary function
summary_func = align_summary_diploid

# Build a dict of stat records

# Collect stats for each region
if aln_input_exists:
    # Alignment file exists, get stats

    stat_list = [
        summary_func(get_len_list(window, input.aln, "subseqfa"))
        for window in region_bed["WINDOW"]
    ]

else:
    # Alignment file does not exist, get empty lists

    stat_list = [summary_func(list())] * region_bed.shape[0]

# Merge stats
if len(stat_list) > 0:
    # Merge summary records
    df = pd.DataFrame(pd.concat(stat_list, axis=1)).T.reset_index(drop=True)

else:
    # Create an empty DataFrame with the right headers if there were no variants in the input
    df = pd.DataFrame([], columns=summary_func(list()).index)

df["HAS_ALN"] = aln_input_exists

# Convert frame types
dtype_dict = {
    field: dtype
    for field, dtype in ALIGN_SUMMARY_FIELD_DTYPE.items()
    if field in df.columns
}

df = df.astype(dtype_dict)

# Prepend region information
df = pd.concat([region_bed, df], axis=1)

# Write
df.to_csv(snakemake.output.tsv, sep="\t", index=False, compression="gzip")
