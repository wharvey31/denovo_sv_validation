


rule callable_bed:
    input:
        bed=find_callable,
        tab=find_bed,
    output:
        tab="temp/CALLABLE/{val_type}/{vartype}_{svtype}/{sample}_{parent}_{hap}.tab",
    resources:
        mem=4,
        hrs=24,
    threads: 1
    script:
        "../scripts/callable_bed.py"


rule combine_callable:
    input:
        all_haps=gather_callable_haps,
    output:
        raw="temp/validation/CALLABLE/{val_type}/{vartype}_{svtype}/{sample}_raw.tsv",
    resources:
        mem=4,
        hrs=24,
    threads: 1
    run:
        for i, file in enumerate(input.all_haps):
            if i == 0:
                df = pd.read_csv(file, sep="\t")
            else:
                df = df.merge(pd.read_csv(file, sep="\t"))
        df.to_csv(output.raw, sep="\t", index=False)


# call_bed_h1 = '/net/eichler/vol27/projects/hprc/nobackups/svpop/hgsvc_hprc_merge/hg38/link_data/pav_raw/{sample}/callable/callable_regions_h1_500.bed.gz',
#       call_bed_h2 = '/net/eichler/vol27/projects/hprc/nobackups/svpop/hgsvc_hprc_merge/hg38/link_data/pav_raw/{sample}/callable/callable_regions_h2_500.bed.gz'
