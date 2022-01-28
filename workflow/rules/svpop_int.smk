


rule int_check:
	input:
		bed = find_bed,
		int_val = find_int
	output:
		'temp/validation/SVPOP/{val_type}/{vartype}_{svtype}/{sample}_raw.tsv'
	run:
		val_df = pd.read_csv(input.int_val, sep='\t')
		all_df = pd.read_csv()