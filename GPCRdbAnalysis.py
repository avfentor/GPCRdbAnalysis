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
    color_plot_no = coloring_plot(full_table)
    saved_graph = draw_graph(full_table, color_plot_no)
    datatable = create_table_for_analysis(color_plot_no)

    return full_table, saved_graph, datatable


def coloring_plot(full_table):
    '''
    Coloring_plot takes the baove created full_table which is unique to each protein family. Based on the
    parameters given in the search, those individual residues will be colored on the graph that is created below.
    A true/false column is being added to the already existing dataframe. We can either give it
    individual residues or parameters.
    '''

    search = choosing_umd_residues(full_table)
    in_search = full_table['new_label'].isin(search)
    TF_table = pd.DataFrame(in_search)
    color_plot = full_table.join(TF_table, lsuffix='_x', rsuffix='_y', )
    color_plot.rename(columns={'new_label_x': 'new_label', 'new_label_y': 'T_or_F'}, inplace=True)
    color_plot['id_if_selected'] = ''
    color_plot.loc[color_plot['T_or_F'], 'id_if_selected'] = color_plot[color_plot['T_or_F']]['new_label']

    occremove = color_plot[color_plot.Occupancy < color_plot.Occupancy.quantile(.20)]['Occupancy'].unique()
    color_plot_no = color_plot[~color_plot['Occupancy'].isin(occremove)]

    return color_plot_no


def draw_graph(full_table, color_plot_no, fix_axes=False):
    '''
    Creates a data table that the graph will be created from. At the end, the graph
    is being saved to the same folder the dataframe came from.
    '''

    # Creating the data table for the plot
    data_plot = pd.melt(color_plot_no,
                        id_vars=color_plot_no.columns[2:].tolist(),
                        var_name='Variant_Effect', value_name='Count')

    # Plotting the above data
    lm = sns.lmplot(x='shenkin', y='Count', col='Variant_Effect', hue="id_if_selected",
                    palette='bright', data=data_plot,
                    fit_reg=False, sharex=True, sharey=True)
    axes = lm.axes
    if fix_axes:
        axes[0, 0].set_xlim(0, fix_axes[0])
        axes[0, 1].set_xlim(0, fix_axes[0])
        axes[0, 0].set_ylim(0, fix_axes[1])
        axes[0, 1].set_ylim(0, fix_axes[1])

    # Adds labels to the colored residues
    data_plot.loc[:, 'Count'] = data_plot.loc[:, 'Count'].fillna(0)
    for p in zip(data_plot['shenkin'], data_plot['Count'], data_plot['id_if_selected'], data_plot['Variant_Effect']):
        ax = axes[0, 1]
        if p[3] == 'missense_variant':
            ax = axes[0, 0]
        ax.text(p[0], p[1], p[2])

    # saves the graph to the respective folder
    plt.suptitle(path + ' graph', x=0.27)  # giving the graph title
    plt.subplots_adjust(top=0.88)  # ensures that the title is saved with the graph
    saved_graph = plt.savefig('./data/' + path + '/colored_graph.png')
    return saved_graph


def create_table_for_analysis(color_plot_no):
    '''
    Creates data table that will be used for analysis. Contains column number that will be used
    in Jalview. IMPORTANT!!! 1 needs to be added to each column number when picking residues
    in Jalview because Python starts counting from 0 and Jalview starts from 1.
    data table is being saved to the folder the data came from.
    '''
    true_values = color_plot_no[color_plot_no['T_or_F']]['new_label']
    color_plot_no.loc[color_plot_no['T_or_F'], 'id_if_selected'] = true_values #selected new_label columns based on T/F conditions
    true_table = pd.DataFrame(true_values)
    true_table.reset_index(inplace=True) #creating new index, so columns are part of the saved table
    selected_columns = true_table.merge(color_plot_no.loc[:,['new_label',
                                                          'missense_variant', 'shenkin', 'Occupancy']],
                                        how='left',on='new_label') #merging desired columns
    newcol = ''
    datatable = selected_columns.assign(Effect=newcol)
    datatable.to_csv('./data/' + path + '/' + path + 'colres_datatable.csv', sep=',', mode='w')
    return datatable



if __name__ == '__main__':
