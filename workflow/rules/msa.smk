manifest_df = pd.read_csv(
    "/net/eichler/vol27/projects/autism_genome_assembly/nobackups/analysis/denovo_validations/sv/denovo_validations_14455.p1_sv_ins.tsv",
    sep="\t",
    index_col="ID",
)

manifest_df["REGION"] = (
    manifest_df["#CHROM"]
    + ":"
    + (manifest_df["POS"] - 1000).astype(str)
    + "-"
    + (manifest_df["END"] + 1000).astype(str)
)


wildcard_constraints:
    hap="hap1|hap2",
    ids="|".join(manifest_df.index),
    sample="14455.p1|14455.fa|14455.mo",


rule all:
    input:
        expand(
            "multialign/results/14455.p1_{ids}_{hap}.shared.bed",
            ids=manifest_df.index,
            child="14455.p1",
            hap=["hap1", "hap2"],
        ),


rule extract_seq:
    input:
        bam="{sample}_{hap}_sorted.bam",
    output:
        fa="multialign/{ids}/{sample}_{hap}.out",
    params:
        region=find_region,
    shell:
        """
        module load subseqfa/1.0
        subseqfa -r {params.region} {input.bam} > {output.fa}
        """


rule rename:
    input:
        fa=rules.extract_seq.output.fa,
    output:
        clean="multialign/{ids}/{sample}_{hap}.out.fa",
    run:
        with open(output.clean, "w") as outfile:
            with open(input.fa, "r") as infile:
                for line in infile:
                    if line[0] == ">":
                        outfile.write(f">{wildcards.sample}_{wildcards.hap}\n")
                    else:
                        outfile.write(line)


rule clustalo:
    input:
        fa=combine_fasta,
    output:
        clust="multialign/results/{ids}/clustal.out",
    shell:
        """
        module load clustal/1.2.4
        cat {input.fa} | clustalo -i - -o {output.clust}
        """


rule process_align:
    input:
        clust=rules.clustalo.output.clust,
    output:
        bed="multialign/results/{ids}/{sample}_{hap}_gap.bed",
    run:
        record_dict = {}
        with open(input.clust) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                text = str(record.seq)
                hap = record.id
                if hap == f"{wildcards.sample}_{wildcards.hap}":
                    for m in re.finditer("-", text):
                        if hap not in record_dict:
                            record_dict[hap] = pd.DataFrame.from_dict(
                                {
                                    "chr": [wildcards.ids],
                                    "start": [m.start()],
                                    "end": [m.end()],
                                    "hap": [hap],
                                }
                            )
                        else:
                            record_dict[hap] = record_dict[hap].append(
                                pd.DataFrame.from_dict(
                                    {
                                        "chr": [wildcards.ids],
                                        "start": [m.start()],
                                        "end": [m.end()],
                                        "hap": [hap],
                                    }
                                )
                            )

        out_df = pd.DataFrame()
        for record in record_dict:
            bed_df = (
                BedTool.from_dataframe(
                    record_dict[record].sort_values(["chr", "start"])
                )
                .merge()
                .to_dataframe(header=None, names=["#chrom", "start", "end"])
            )
            bed_df["hap"] = record
            out_df = out_df.append(bed_df)
        out_df.to_csv(output.bed, sep="\t", index=False)


rule check_hap:
    input:
        hap="multialign/results/{ids}/14455.p1_{hap}_gap.bed",
        parents=expand(
            "multialign/results/{{ids}}/{sample}_{hap}_gap.bed",
            hap=["hap1", "hap2"],
            sample=["14455.fa", "14455.mo"],
        ),
    output:
        bed="multialign/results/14455.p1_{ids}_{hap}.shared.bed",
    run:
        try:
            df = pd.read_csv(input.hap, sep="\t")
        except:
            df = pd.DataFrame(columns=["#CHROM", "POS", "END", "sample"])

        a = BedTool(input.parents[0])
        b = BedTool(input.parents[1])
        c = BedTool(input.parents[2])
        d = BedTool(input.parents[3])

        par_bed = a.intersect(b).intersect(c).intersect(d)

        shared_ovl = (
            BedTool.from_dataframe(df)
            .intersect(par_bed)
            .to_dataframe(header=None, names=["#CHROM", "POS", "END", "sample"])
        )

        if len(shared_ovl) > 0:
            shared_ovl["ovl"] = shared_ovl["END"] - shared_ovl["POS"]

        shared_ovl.to_csv(output.bed, sep="\t", index=False)


# rule overlap_self:
#     input:
#         h1 = expand('multialign/results/14455.p1_{{ids}}_{hap}.shared.bed', hap=['hap1']),
#         h2 = expand('multialign/results/14455.p1_{{ids}}_{hap}.shared.bed', hap=['hap2']),
#     output:
#         val = '14455.p1_val/{ids}.val'
#     run:
#         if len(BedTool.from_dataframe(pd.read_csv(input.h1, sep='\t')))
#         a = BedTool(input.bed)
#         # b = BedTool(input.bed)
#         df = a.intersect(b, wb=True).to_dataframe(header=None,names=['#ID', 'start_1', 'end_1', 'name_1', 'chr_2', 'start_2', 'end_2', 'name_2'])
#         # df = df.loc[(df['name_1'].str.contains('14455.p1')) & (~df['name_2'].str.contains('14455.p1'))].groupby(['name_1', 'name_2', 'chr_1']).count().reset_index()
#         df = df.loc[(df['name_1'].str.contains('14455.p1')) & (~df['name_2'].str.contains('14455.p1'))]
#         df['ovl'] = df['end_1'].astype(int) - df['start_1'].astype(int)
#         df[['ID', 'name_1', 'name_2', 'ovl', 'start_1', 'end_1']].to_csv(output.val, sep='\t', index=False)
#         # df
