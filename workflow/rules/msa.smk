


rule extract_seq:
    input:
        bam=find_asm_aln,
    output:
        fa="temp/validation/ASM/{val_type}/{sample}/{vartype}_{svtype}/{ids}_{hap}.out",
    params:
        region=find_region,
    resources:
        mem=4,
        hrs=24,
    threads: 1
    shell:
        """
        subseqfa -r {params.region} {input.bam} > {output.fa}
        """


rule rename:
    input:
        fa=rules.extract_seq.output.fa,
    output:
        clean="temp/validation/ASM/{val_type}/{sample}/{vartype}_{svtype}/{ids}_{hap}.out.fa",
    resources:
        mem=4,
        hrs=24,
    threads: 1
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
        clust="temp/validation/ASM/{val_type}/{sample}/clustalo/{vartype}_{svtype}/{ids}/clustal.out",
    resources:
        mem=lambda wildcards, attempt : 4**attempt,
        hrs=24,
    threads: 1
    shell:
        """
        if [[ $( grep ">" {input.fa} | wc -l ) < 2 ]]; then
            touch {output.clust}
        else
            cat {input.fa} | clustalo -i - -o {output.clust}
        fi
        """


rule process_align:
    input:
        clust=rules.clustalo.output.clust,
    output:
        bed="temp/validation/ASM/{val_type}/{vartype}_{svtype}/{ids}/{sample}_{hap}_gap.bed",
    resources:
        mem=4,
        hrs=24,
    threads: 1
    run:
        record_dict = {}
        with open(input.clust) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                text = str(record.seq)
                hap = record.id
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
        hap="temp/validation/ASM/{val_type}/{vartype}_{svtype}/{ids}/{sample}_{hap}_gap.bed",
    output:
        bed="temp/validation/ASM/{val_type}/{vartype}_{svtype}/{ids}/{sample}_{hap}_shared.bed",
    resources:
        mem=4,
        hrs=24,
    threads: 1
    run:
        try:
            df = pd.read_csv(input.hap, sep="\t")
        except:
            df = pd.DataFrame(columns=["#chrom", "start", "end", "hap"])

        df_hap = df.loc[df["hap"] == f"{wildcards.sample}_{wildcards.hap}"].copy()

        df_parents = df.loc[df["hap"] != f"{wildcards.sample}_{wildcards.hap}"].copy()

        par_bed = BedTool.from_dataframe(df_parents)

        for i, sample in enumerate(df_parents["hap"].unique()):
            if i == 0:
                par_bed = BedTool.from_dataframe(
                    df_parents.loc[df_parents["hap"] == sample]
                )
            else:
                par_bed = par_bed.intersect(
                    BedTool.from_dataframe(df_parents.loc[df_parents["hap"] == sample])
                )

        shared_ovl = (
            BedTool.from_dataframe(df)
            .intersect(par_bed)
            .to_dataframe(header=None, names=["#CHROM", "POS", "END", "sample"])
        )

        if len(shared_ovl) > 0:
            shared_ovl["ovl"] = shared_ovl["END"] - shared_ovl["POS"]

        shared_ovl.to_csv(output.bed, sep="\t", index=False)


rule combine_ids:
    input:
        bed=find_ids,
        all_bed=find_bed,
    output:
        val="temp/validation/ASM/{val_type}/{vartype}_{svtype}/{sample}_raw.tsv",
    resources:
        mem=16,
        hrs=24,
    threads: 1
    run:
        out_df = pd.DataFrame()
        df = pd.concat(
            [
                pd.read_csv(
                    file,
                    sep="\t",
                    header=0,
                    names=["#CHROM", "POS", "END", "sample", "ovl"],
                )
                for file in input.bed
            ]
        )
        df = df.drop_duplicates().dropna()
        df["SVLEN"] = df["#CHROM"].str.split("-", expand=True)[3].astype(int)
        df["LEN_MATCH"] = df.apply(
            lambda row: "YES" if row["SVLEN"] == row["ovl"] else "NO", axis=1
        )
        df = df.loc[df["LEN_MATCH"] == "YES"]
        for svid in df["#CHROM"].unique():
            sv_df = df.loc[df["#CHROM"] == svid]
            dn_sample = [
                x
                for x in [f"{wildcards.sample}_hap1", f"{wildcards.sample}_hap2"]
                if x not in sv_df["sample"].values
            ]
            if len(dn_sample) == 1:
                val_call = "VALID"
                val_sample = dn_sample[0]
            else:
                val_call = "NOTVALID"
                val_sample = ""
            out_df = out_df.append(
                pd.DataFrame.from_dict(
                    {
                        f"MSA_VAL_SAMPLE_{wildcards.val_type}": [val_sample],
                        f"MSA_VAL_{wildcards.val_type}": [val_call],
                        "ID": [svid],
                    }
                )
            )

        bed_df = pd.read_csv(input.all_bed, sep="\t", usecols=["ID"])

        out_df = out_df.merge(bed_df, how="outer")

        out_df[
            [
                "ID",
                f"MSA_VAL_{wildcards.val_type}",
            ]
        ].to_csv(output.val, sep="\t", index=False)
