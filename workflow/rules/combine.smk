






rule minimize_input:
	input:
		file = 'temp/validation/{val_type}/{val_method}/{sample}_raw.tsv'
	output:
		val = 'temp/validation/{val_type}/{val_method}/{sample}_val.tsv'
	threads: 1
	resources:
		mem = 16,
		hrs = 12
	run:
		if val_method == 'SVPOP':
			
		elif val_method == 'CALLABLE':

		else:

			







rule find_combine:
	input:
		val = determine_input
	output:
		val = 'dn_val/{sample}_valid_table.tsv',
		all_calls = 'dn_val/{sample}_all_table.tsv'
	threads: 1
	resources:
		mem = 16,
		hrs = 12
	run:
		for i, file in enumerate(input.val):
			if i == 0:
				df = pd.read_csv(file, sep='\t')
			else:
				df = df.merge(pd.read_csv(file, sep='\t'))
		sample_cols = [x for x in df.columns if wildcards.sample in x and x != 'ID']
		parental_cols = [x for x in df.columns if wildcards.sample not in x and x != 'ID']
		sample_val = df[sample_cols].str.replace('VALID', True)
		parent_val = df[parental_cols].str.replace('NOTVALID', True)
		df.to_csv(output.all_calls, sep='\t', index=False)
		sample_val['ID'].merge(parent_val['ID']).to_csv(output.val, sep='\t', index=False)




