


rule callable_bed:
    input:
        bed=find_callable,
        tab=find_calls,
    output:
        tab="temp/CALLABLE/{parent}_{hap}.tab",
    run:
        df = pd.read_csv(input.tab, sep="\t")
        call = Bedtool(input.bed)
        regions = (
            BedTool.from_dataframe(df[["#CHROM", "POS", "END", "ID"]])
            .intersect(call, wa=True)
            .to_dataframe(names=["#CHROM", "POS", "END", "ID"])
        )
        regions[f"CALLABLE_{wildcards.parent}_{wildcards.hap}"] = "VALID"
        df = df.merge(regions, how="left")
        df[f"CALLABLE_{wildcards.parent}_{wildcards.hap}"] = df[
            f"CALLABLE_{wildcards.parent}_{wildcards.hap}"
        ].fillna("NOTVALID")
        df.to_csv(output.tab, sep="\t", index=False)


rule combine_callable:
	input:
        


# call_bed_h1 = '/net/eichler/vol27/projects/hprc/nobackups/svpop/hgsvc_hprc_merge/hg38/link_data/pav_raw/{sample}/callable/callable_regions_h1_500.bed.gz',
# 		call_bed_h2 = '/net/eichler/vol27/projects/hprc/nobackups/svpop/hgsvc_hprc_merge/hg38/link_data/pav_raw/{sample}/callable/callable_regions_h2_500.bed.gz'
