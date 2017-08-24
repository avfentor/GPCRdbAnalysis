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


def _parse_display_generic_number(x):
    """
    Parse "display_generic_number".

    GPCRdb numbering looks like "1.26x26". This function splits each component out of a list
    of such numbers and returns a DataFrame.
    """
    parsed = pd.Series(x).str.split('[.x]', expand=True)
    parsed.columns = ['seqment', 'conventional', 'GPCRdb']

    return parsed


def fetch_gpcr_db_lookup_table(protein_name, parse_dgn=True):
    """
    Retrieve GPCRdb residue numbering lookup table.
    """
    # GET request
    url = 'http://gpcrdb.org/services/residues/{0}/'.format(protein_name)
    response = requests.get(url)

    # Format response into dataframe
    protein_mapping = pd.DataFrame(response.json())
    protein_mapping.insert(0, 'SOURCE_ID', protein_name)

    # Parse diplay_generic_number
    if parse_dgn:
        reformed = _parse_display_generic_number(protein_mapping['display_generic_number'])
        protein_mapping = pd.concat(axis=1, objs=[protein_mapping, reformed])

    return protein_mapping


def filter_empty_sequences(alignment):
    #Check for empty sequences
    passed = []
    for seq in alignment:
        sequence_string = str(seq.seq)
        if not all([x == '-' for x in sequence_string]):
            passed.append(seq)

    new_alignment = MultipleSeqAlignment(passed)
    return new_alignment


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


    # Plotting results
    saved_graph = draw_graph(full_table, color_plot)

    return full_table, saved_graph


def draw_graph(full_table, color_plot, fix_axes=False):
    '''
    Creates a data table that the graph will be created from. At the end, the graph
    is being saved to the same folder the dataframe came from.
    '''

    # Creating the data table for the plot
    data_plot = pd.melt(color_plot,
                        id_vars=color_plot.columns[2:].tolist(),
                        var_name='Variant_Effect', value_name='Count')

    # Plotting the above data
    lm = sns.lmplot(x='shenkin', y='Count', col='Variant_Effect', hue="Variant_Effect",
                    palette='bright', data=data_plot,
                    fit_reg=True, sharex=True, sharey=True)
    axes = lm.axes
    if fix_axes:
        axes[0, 0].set_xlim(0, fix_axes[0])
        axes[0, 1].set_xlim(0, fix_axes[0])
        axes[0, 0].set_ylim(0, fix_axes[1])
        axes[0, 1].set_ylim(0, fix_axes[1])

    # saves the graph to the respective folder
    plt.suptitle(path + ' graph', x=0.27)  # giving the graph title
    plt.subplots_adjust(top=0.88)  # ensures that the title is saved with the graph
    saved_graph = plt.savefig('./data/' + path + '/_graph.png')
    return saved_graph



if __name__ == '__main__':
