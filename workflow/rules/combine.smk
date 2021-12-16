



rule minimize_input:
	input:
		file = 'temp/validation/{val_method}/{val_type}/{sample}_raw.tsv'
	output:
		val = 'temp/validation/{val_method}/{val_type}/{sample}_val.tsv'
	threads: 1
	resources:
		mem = 16,
		hrs = 12
	run:
		base_columns = ['']
	# 	if val_type == 'SVPOP':

		elif val_method == 'CALLABLE':


	# 	else:


rule find_combine:
    input:
        val=determine_input,
    output:
        val="dn_val/{sample}_valid_table.tsv",
        all_calls="dn_val/{sample}_all_table.tsv",
    threads: 1
    resources:
        mem=16,
        hrs=12,
    script:
        "../scripts/combine.py"


