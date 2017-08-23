import argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment #reading in the alignment
from scipy.stats import linregress
from varalign import jalview as jv
import pandas as pd #dataframe
import os
import requests #getting GPCR lookup table from GPCRdb website
import requests_cache
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt

%matplotlib inline

requests_cache.install_cache('gpcrdb_cache')


def calculate_occupancy(new_alignment):
    #Calculating occupancy in alignment
    gaps = []
    occupancy = []
    alignment_count = len(new_alignment)
    alignment_length = new_alignment.get_alignment_length()
    for a in range(alignment_length):
            column_count = new_alignment[:, a]
            gap_count = column_count.count('-')
            gaps.append(gap_count)
            occupancy.append(alignment_count - gap_count)
    occupancy_table = pd.DataFrame(occupancy)
    occupancy_table.index.name = 'column'
    occupancy_table.columns = ['Occupancy']
    occupancy_table.reset_index(inplace=True)
    return occupancy_table


def choosing_umd_residues(full_table):
    #creating dataframe using only the values we need (missense and synonymous variant, shenkin score occupancy, new_label)
    scores = full_table[full_table.columns[0:2]]
    score = full_table[full_table.columns[5:6]]
    info = full_table[full_table.columns[20:22]]
    required_data = scores.join([score,info], how='outer')

    #Ranking data from high to low based on shenkin score
    sorted_data = required_data.sort_values(['shenkin', 'missense_variant'],
                                            axis=0, ascending=False, inplace=False,
                                            kind='quicksort', na_position='last')
    #Choosing only the top 95% of shenkin scores
    shenkin_data = sorted_data[sorted_data.shenkin > sorted_data.shenkin.quantile(.95)]
    #Filtering out high missense variant counts
    umd_residues = shenkin_data[shenkin_data.missense_variant < shenkin_data.missense_variant.quantile(.90)]['new_label'].unique()
    return umd_residues


def pipeline(path):
    """
    The below code imports the fasta file and csv file from the folders given above,
    checks for gaps and occupancy,and merges variants on table.
    Gives each residue a new label based on their TM position and residue number.
    Returns new_alignment, conservation_plane and alignment_columns (alignment, 2 tables),
    which are used later in the analysis.
    """
    # Reading the CSV and FASTA file from given location
    GPCR_alignment_csv = os.path.join(path, 'GPCRdb_alignment.csv')
    GPCR_alignment_fasta = os.path.join(path, 'GPCRdb_alignment.fasta')

    # Define AACon output file name
    GPCR_alignment_aacon = GPCR_alignment_fasta.replace('fasta', 'aacons')
    GPCR_alignment_aacon = GPCR_alignment_aacon.replace('fa', 'aacons')

    # Fixed FASTA file
    GPCR_alignment_fixed = GPCR_alignment_fasta.replace('.', '_fixed.')

    print
    'Reading {} file'.format(path)

    # Reading CSV table with Pandas
    csv_table = pd.read_csv(GPCR_alignment_csv, header=1, index_col=0)
    csv_table.reset_index(inplace=True)
    seq_map = pd.melt(csv_table, id_vars='index')
    seq_map.columns = ['seq_name', 'gpcr_residue_number', 'aa']

    # Parse alignment column GPCRdb IDs
    alignment_columns = _parse_display_generic_number(csv_table.columns[1:])
    alignment_columns.index.name = 'column'
    alignment_columns.reset_index(inplace=True)
    new_label = alignment_columns['seqment'].str.cat(alignment_columns['GPCRdb'],
                                                     sep='x')
    alignment_columns.insert(1, 'new_label', new_label)

    # Reading the alignment from the alignment_path
    alignment = AlignIO.read(GPCR_alignment_fasta, format='fasta')

    # Checking for empty sequences and creating a fixed alignment file
    new_alignment = filter_empty_sequences(alignment)
    with open(GPCR_alignment_fixed, 'wb') as out:
        AlignIO.write(new_alignment, out, 'fasta')


        # Merge sequences onto alignment columns
    alignment_mapped_proteins = []
    for seq in new_alignment:
        protein_mapping_table = fetch_gpcr_db_lookup_table(seq.name)
        mapped = pd.merge(alignment_columns, protein_mapping_table,
                          how='inner', on=['seqment', 'GPCRdb'])
        alignment_mapped_proteins.append(mapped)

    # Joins all individual sequence mapped tables
    mapped_alignment = pd.concat(alignment_mapped_proteins)

    # Reading the variant table from CSV file
    variant_table = pd.read_csv('gpcr_db_class_a_variants_filtered.csv', index_col=0, low_memory=False)
    variant_table['Protein_position'] = pd.to_numeric(variant_table['Protein_position'], errors='coerce')

    # Merging the alignment table and the variant table together
    alignment_variant_table = pd.merge(mapped_alignment, variant_table,
                                       left_on=['SOURCE_ID', 'sequence_number'],
                                       right_on=['SOURCE_ID', 'Protein_position'])
    # Variant count
    conseq_groups = alignment_variant_table.groupby(['column',
                                                     'Consequence']).size().unstack(fill_value=0)
    # Reading FASTA alignment and creating AACONS file
    subprocess.call(["java", "-jar", "compbio-conservation-1.1.jar",
                     "-i=" + GPCR_alignment_fixed,
                     "-o=" + GPCR_alignment_aacon])

    # Reading AACons table with pandas
    aacons_table = pd.read_csv(GPCR_alignment_aacon, sep=' ',
                               header=None,
                               index_col=0).transpose().dropna()
    aacons_table.columns = aacons_table.columns.str.replace('#', '')
    aacons_table.columns = aacons_table.columns.str.lower()
    aacons_table.index.name = 'column'
    aacons_table.index = aacons_table.index - 1

    # AACons and column stat joint table
    conservation_plane = conseq_groups[['missense_variant', 'synonymous_variant']].join(aacons_table, how='outer')

    # Calculating occupancy and making full table
    occupancy_table = calculate_occupancy(new_alignment)
    scores = conservation_plane.join(occupancy_table).drop(['column'], axis=1)
    full_table = scores.join(alignment_columns).drop(['column'], axis=1)

    # Plotting results
    color_plot = coloring_plot(full_table)
    saved_graph = draw_graph(full_table, color_plot)
    datatable = create_table_for_analysis(color_plot)

    return full_table, saved_graph, datatable


