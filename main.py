import numpy as np
from csv import reader as csv_reader
from matplotlib import pyplot as plt
import logging
import sys
from scipy.stats import gaussian_kde, norm, chi2
import pickle
from itertools import product
from collections import defaultdict


def set_up_loggers():
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    root.addHandler(ch)

    ch1 = logging.FileHandler('info.log')
    ch1.setLevel(logging.INFO)
    ch1.setFormatter(formatter)
    root.addHandler(ch1)

    ch2 = logging.FileHandler('debug.log')
    ch2.setLevel(logging.DEBUG)
    ch2.setFormatter(formatter)
    root.addHandler(ch2)


# define a helper function that manages nans
def nan_helper(value):
    try:
        output = float(value)
    except Exception as error:
        logging.debug(error)
        output = np.nan
    return output


# perform the actual main reading loop
def parse_weights(weights_data_file, transcription_factors_line, start_of_table_line, end_of_table):
    deleted_genes = []
    read_out_genes = []
    _data_matrix = []

    with open(weights_data_file, 'rb') as in_file:
        reader = csv_reader(in_file, delimiter='\t')
        for line_number, line in enumerate(reader):
            if line_number == transcription_factors_line:
                deleted_genes = [string.split(' ')[0] for string in line[1:]]
            if end_of_table > line_number > start_of_table_line:
                read_out_genes.append(line[0])
                # print line
                # print line_number
                _data_matrix.append([nan_helper(value) for value in line[1:]])
                # one line is all experimental conditions for a single gene

    _data_matrix = np.array(_data_matrix)

    # check that we got everything properly
    logging.info('deleted genes:\t%s', len(deleted_genes))
    logging.info('genes in read-out:\t%s', len(read_out_genes))
    logging.info('data matrix shape:\t%s', _data_matrix.shape)



    return _data_matrix, deleted_genes, read_out_genes


def parse_direct_connections(direct_connections_data_file):
    TF_2_Genes = defaultdict(lambda: 0)
    with open(direct_connections_data_file) as source:
        reader = csv_reader(source, delimiter=';')
        for line in reader:
            TF_2_Genes[line[0], line[1]] = 1

    return TF_2_Genes


def parse_names_maps(names_maps_files):
    ens_2_names = {}
    names_2_ens = {}
    with open(names_maps_files) as source:
        reader = csv_reader(source, delimiter="\t")
        for line in reader:
            print line, len(line)
            if line[1]:
                ens_2_names[str.upper(line[0])] = str.upper(line[1])
                names_2_ens[str.upper(line[1])] = str.upper(line[0])

    return ens_2_names, names_2_ens


def calculate_raw_activations(names_map, reverse_names_map, structure, weights, TFs, genes):
    TFs_names = []
    TF_trans_counter = [0, 0]

    genes_names = []
    genes_trans_counter = [0, 0]

    for TF in TFs:
        TF_name = TF.split(' ')[0]
        try:
            TFs_names.append(names_map[str.upper(TF_name)])
        except KeyError:
            TF_trans_counter[1] += 1
            # print 'Failed to translate TF: %s' % TF_name
        else:
            TF_trans_counter[0] += 1
            # print 'Translated TF: %s to %s' % (TF_name, names_map[str.upper(TF_name)])

    for gene in genes:
        try:
            genes_names.append(names_map[str.upper(gene)])
            # print 'Translated gene: %s to %s' % (gene, names_map[str.upper(gene)])
        except KeyError:
            genes_trans_counter[1] += 1
            # print 'Failed to translate gene: %s' % gene
        else:
            genes_trans_counter[0] += 1
            # print 'Translated gene: %s to %s' % (gene, names_map[str.upper(gene)])

    print "TFs translation: %s successes; %s failures" % tuple(TF_trans_counter)
    print "Genes translation: %s successes; %s failures" % tuple(genes_trans_counter)

    TFs_names = []
    TF_trans_counter = [0, 0]

    genes_names = []
    genes_trans_counter = [0, 0]

    separate_TFs = set()
    separate_genes = set()

    re_structure = defaultdict(lambda: 0)

    for TF, gene in structure:
        separate_TFs.add(TF)
        separate_genes.add(gene)
        try:
            re_structure[(reverse_names_map[str.upper(TF)],
                          reverse_names_map[str.upper(gene)])]
        except KeyError:
            pass


    for TF in separate_TFs:
        try:
            TFs_names.append(reverse_names_map[str.upper(TF)])
        except KeyError:
            TF_trans_counter[1] += 1
            # print 'Failed to translate TF: %s' % TF
        else:
            TF_trans_counter[0] += 1
            # print 'Translated TF: %s to %s' % (TF, reverse_names_map[str.upper(TF)])

    for gene in separate_genes:
        try:
            genes_names.append(reverse_names_map[str.upper(gene)])
        except KeyError:
            genes_trans_counter[1] += 1
            # print 'Failed to translate gene: %s' % gene
        else:
            genes_trans_counter[0] += 1
            # print 'Translated gene: %s to %s' % (gene, reverse_names_map[str.upper(gene)])

    print "TFs reverse translation: %s successes; %s failures" % tuple(TF_trans_counter)
    print "Genes reverse translation: %s successes; %s failures" % tuple(genes_trans_counter)

    for TF, gene in re_structure:
        if gene in genes and TF is TFs:
            re_structure[(TF, gene)] = weights[TFs.index(TF), genes.index(gene)]

    print len(re_structure)

    return re_structure


if __name__ == "__main__":

    set_up_loggers()

    connections_matrix = parse_direct_connections("RegulationTwoColumnTable_Documented_2013927.tsv")

    names_map, reverse_names_map = parse_names_maps("Ensembl_to_gene_names_map.txt")

    weights_matrix, TFs, genes = parse_weights(
        weights_data_file='GSE4654_series_matrix.txt',  # data source for parsing
        transcription_factors_line=44,
        start_of_table_line=79,
        end_of_table=6509)

    raw_activations = calculate_raw_activations(names_map,
                                                reverse_names_map,
                                                connections_matrix,
                                                weights_matrix,
                                                TFs, genes)

    raise Exception('Debugging')

    # Show overall factor distribution
    data = weights_matrix.flatten()
    fltr = np.logical_not(np.isnan(data))
    density = gaussian_kde(data[fltr].flatten())
    xs = np.linspace(data[fltr].min(), data[fltr].max(), 200)
    plt.title('distribution of gene expression alteration in response to TF deletion')
    ax = plt.plot(xs, -np.sqrt(-np.log(density(xs))), 'k')
    plt.plot(xs, -np.sqrt(-np.log(norm.pdf(xs))), 'r')
    plt.xlabel('z-score deviation')
    plt.ylabel('distribution density (log)')
    plt.show()

    data = np.power(weights_matrix, 2)
    fltr = np.logical_not(np.isnan(data))
    density = gaussian_kde(data[fltr].flatten())
    xs = np.linspace(np.sqrt(data)[fltr].min(), np.sqrt(data)[fltr].max(), 200)
    plt.title('distribution of gene expression alteration in response to TF deletion')
    ax = plt.plot(xs, density(np.power(xs, 2)), 'k')
    plt.semilogy(xs, chi2.pdf(np.power(xs, 2), 1), 'r')
    plt.xlabel('z-score deviation')
    plt.ylabel('distribution density (log)')
    plt.show()

    # Generate a plot for deviation in each experimental condition
    data2 = np.sqrt(np.nanmean(np.power(weights_matrix, 2), axis=1))
    ## per experimental condition
    fltr = np.logical_not(np.isnan(data2))
    density = gaussian_kde(data2[fltr].flatten())
    xs = np.linspace(data2[fltr].min(), data2[fltr].max(), 200)
    plt.title('distribution of euclidean distances of gene expression deviation\nper-condition')
    ax = plt.semilogy(xs, density(xs), 'k')
    plt.xlabel('mean L2 z-score deviation')
    plt.ylabel('distribution density (log)')
    plt.show()

    data3 = np.sqrt(np.nanmean(np.power(weights_matrix, 2), axis=0))
    ## per gene
    fltr = np.logical_not(np.isnan(data3))
    density = gaussian_kde(data3[fltr].flatten())
    xs = np.linspace(data3[fltr].min(), data3[fltr].max(), 200)
    plt.title('distribution of euclidean distances of gene expression deviation\nper-gene')
    ax = plt.semilogy(xs, density(xs), 'k')
    plt.xlabel('mean L2 z-score deviation')
    plt.ylabel('distribution density (log)')
    plt.show()

    # pickle.dump(weights_matrix[np.logical_or(weights_matrix < -2, weights_matrix > 2)], open('activations.dmp', 'wb'))
    pickle.dump(weights_matrix, open('activations.dmp', 'wb'))
