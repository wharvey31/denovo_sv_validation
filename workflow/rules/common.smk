"""
Subseq stats.
"""

# Alignment summary table data types.
ALIGN_SUMMARY_FIELD_DTYPE = {
    'N': np.int32,
    'MEAN': np.float32,
    'MIN': np.int32,
    'MAX': np.int32,
    'N_LO': np.int32,
    'MEAN_LO': np.float32,
    'MIN_LO': np.int32,
    'MAX_LO': np.int32,
    'N_HI': np.int32,
    'MEAN_HI': np.float32,
    'MIN_HI': np.int32,
    'MAX_HI': np.int32,
    'DIST_LH': np.int32,
    'HAS_ALN': bool
}


SET_DEF = {
    'indel20-50':
        (20, 50, 5)
    'sv50-100':
        (50, 100, 20),
    'sv100-200':
        (100, 200, 25),
    'sv200-500':
        (200, 500, 40),
    'sv500-1000':
        (500, 1000, 100),
    'sv1-2k':
        (1000, 2000, 200),
    'sv2-4k':
        (2000, 4000, 250),
    'sv4k-max':
        (4000, None, 300),
}

def find_bed(wildcards):
    return manifest_df.at[wildcards.sample, 'BED']


#
# StepMiner
#

def step_miner(len_list):
    """
    Take a sorted list and split it into two sets (low and high) for each possible split (index = 1, 2, 3, ..., len(len_list)).
    For each set (low and high), compute the root-mean-squared error between the mean of the set and each element. For the split with
    the lowest error, the low and high lists are returned (tuple of two lists, low list is the first element, high list is the second).

    This StepMiner algorithm is borrowed from Debashis Sahoo's thesis:
    http://genedesk.ucsd.edu/home/dsahoo-thesis.pdf
    """

    min_index = None
    min_error = None

    n_list = len(len_list)

    if len(len_list) == 0:
        return [], [], 0, 0

    if len(len_list) == 1:
        return len_list, [], 0, 0

    # Get split with the least error
    for index in range(1, len(len_list)):

        len_low = len_list[:index]
        len_high = len_list[index:]

        mean_low = np.mean(len_low)
        mean_high = np.mean(len_high)

        error = (
            np.sum(np.abs(len_low - mean_low)**2) +
            np.sum(np.abs(len_high - mean_high)**2)
        )

        if error > 0:
            error = np.log2(error / n_list)

        if index == 1 or error < min_error:
            min_error = error
            min_index = index

    # Return splits
    return len_list[:min_index], len_list[min_index:], min_index, min_error

#
# Alignment summary record
#

def align_summary_haploid(len_list):
    
    # Get stats
    n = len(len_list)

    if n > 0:
        mean = np.mean(len_list)
        min = np.min(len_list)
        max = np.max(len_list)
    else:
        mean = 0
        min = 0
        max = 0

    # Return series
    return pd.Series(
        [
            n, mean, min, max,
            ','.join(['{:d}'.format(val) for val in len_list])
        ],
        index=[
            'N', 'MEAN', 'MIN', 'MAX',
            'LENS'
        ]
    )

def align_summary_diploid(len_list):

    len_low, len_high, min_index, min_error = step_miner(sorted(len_list))

    # Get stats
    n_low = len(len_low)
    n_high = len(len_high)

    if n_low > 0:
        mean_low = np.mean(len_low)
        min_low = np.min(len_low)
        max_low = np.max(len_low)
    else:
        mean_low = 0
        min_low = 0
        max_low = 0

    if n_high > 0:
        mean_high = np.mean(len_high)
        min_high = np.min(len_high)
        max_high = np.max(len_high)
    else:
        mean_high = 0
        min_high = 0
        max_high = 0

    if n_low > 0 and n_high > 0:
        separation = min_high - max_low
    else:
        separation = 0

    # Return series
    return pd.Series(
        [
            n_low, mean_low, min_low, max_low,
            n_high, mean_high, min_high, max_high,
            separation,
            ','.join(['{:d}'.format(val) for val in len_low]),
            ','.join(['{:d}'.format(val) for val in len_high])
        ],
        index=[
            'N_LO', 'MEAN_LO', 'MIN_LO', 'MAX_LO',
            'N_HI', 'MEAN_HI', 'MIN_HI', 'MAX_HI',
            'DIST_LH',
            'LENS_LO',
            'LENS_HI'
        ]
    )


"""
Validation functions.
"""

def validate_summary(df, strategy='size50_2_4'):
    """
    Run validation on. Generates a "VAL" colmun with:
    * VALID: 
    * NOTVALID: 
    * NOCALL: 
    * NODATA: 
    """
    
    # Get parameters
    match_obj = re.match(r'^size(\d+)_(\d+)_(\d+)$', strategy)
    
    if match_obj is None:
        raise RuntimeError('No implementation for strategy: ' + strategy)
    
    val_threshold = np.float32(match_obj[1]) / 100
    
    min_support = np.int32(match_obj[2])
    min_call_depth = np.int32(match_obj[3])
    
    # Subset df to needed columns
    df = df[['ID', 'SAMPLE', 'CALLER', 'ALNSAMPLE', 'ALNSOURCE', 'SVTYPE', 'SVLEN', 'LENS_HI', 'LENS_LO', 'HAS_ALN', 'WINDOW_SIZE']].copy()
    
    # Get lengths and SV length differences
    df['LEN'] = df.apply(lambda row:
        (row['LENS_LO'].split(',') if not pd.isnull(row['LENS_LO']) else []) +
        (row['LENS_HI'].split(',') if not pd.isnull(row['LENS_HI']) else []),
        axis=1
    )
    
    df['LEN_DIFF'] = df.apply(lambda row:
        [
            (int(val) - row['WINDOW_SIZE'] - (
                row['SVLEN'] if row['SVTYPE'] == 'INS' else 0)
            )
            for val in row['LEN']
        ] ,
        axis=1
    )
    
    # Count support
    df['SUPPORT_COUNT'] = df.apply(
        lambda row: np.sum(np.abs([
            np.int32(element) / row['SVLEN'] for element in row['LEN_DIFF']
        ]) < val_threshold),
        axis=1
    )
    
    # Call validation status
    df['VAL'] = df.apply(lambda row:
        (
            'VALID' if (row['SUPPORT_COUNT'] > min_support) else 'NOTVALID'
        ) if (
            len(row['LEN']) >= min_call_depth
        ) else 'NOCALL',
        axis=1
    )
        
    df['VAL'] = df.apply(lambda row: row['VAL'] if row['HAS_ALN'] else 'NODATA', axis=1)
    
    # Clean up
    del(df['LENS_HI'])
    del(df['LENS_LO'])
    del(df['LEN'])
    
    df['LEN_DIFF'] = df['LEN_DIFF'].apply(lambda vals: ','.join([str(val) for val in vals]))
    
    return df


def get_ref_fai(fai_file_name):
    """
    Read a reference FAI file as a Series object of sequence lengths keyed by sequence name.
    """
    
    return pd.read_csv(
        fai_file_name,
        sep='\t',
        names=('CHROM', 'LEN', 'START', 'LEN_BP', 'LEN_BYTES'),
        usecols=('CHROM', 'LEN'),
        index_col='CHROM',
        squeeze=True
    )

def get_lengths(fa_file):
    """
    Get a list of record lengths for one input file.
    """
    with gzip.open(fa_file, 'rt') as in_file:
        return [len(record.seq) for record in SeqIO.parse(in_file, 'fasta')]

def get_len_list(window, aln_input, subseq_exe):
    """
    Get a list of alignment record lengths over a window.
    
    :param window: Position string (chrom:pos-end). Coordinates are 1-based inclusive (not BED).
    :param aln_input: Alignment input as BAM or CRAM.
    :param subseq_exe: Path to subseq executable.
    """

    # Open process
    proc = subprocess.Popen([subseq_exe, '-b', '-r', window, aln_input], stdout=subprocess.PIPE)

    stdout, stderr = proc.communicate()
    
    # Return list
    with io.StringIO(stdout.decode()) as in_file:
        return [len(record.seq) for record in SeqIO.parse(in_file, 'fasta')]

def determine_combined_set(wildcards):
    val_types = [x for x in [config.get('ASM'), config.get('READS'), config.get('SVPOP'), config.get('CALLABLE')] if x != None]
    full_set = []
    for check in val_types:
        for key in val_types.keys():
            full_set.append('temp/validation/{check}/{key}/{sample}_val.tsv'.format(check=check, key=key, sample='{sample}'))
    return full_set



