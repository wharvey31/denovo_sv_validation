

rule combine_val:
    input:
        val=determine_input,
        bed=find_bed,
    output:
        val="dn_val/{sample}_valid_table.tsv",
        all_calls="dn_val/{sample}_all_table.tsv",
    threads: 1
    resources:
        mem=16,
        hrs=12,
    script:
        "../scripts/combine.py"
