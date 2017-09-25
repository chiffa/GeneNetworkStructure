import numpy as np
from csv import reader as csv_reader
from matplotlib import pyplot as plt
import logging
import sys
from scipy.stats import gaussian_kde, norm, chi2
import pickle

# Set up the logger
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

# set up the input and lines of interest anchors
input_file = 'GSE4654_series_matrix.txt'
genes_line = 44
start_of_table_line = 79
end_of_table = 6509

# set up the structures
deleted_genes = []
read_out_genes = []
data_matrix = []


# define a helper function that manages nans
def helper_function(string):
    try:
        output = float(value)
    except Exception as error:
        logging.debug(error)
        output = np.nan
    return output

# perform the actual main reading loop
with open(input_file, 'rb') as in_file:
    reader = csv_reader(in_file, delimiter='\t')
    for line_number, line in enumerate(reader):
        if line_number == genes_line:
            deleted_genes = [string.split(' ')[0] for string in line[1:]]
        if line_number > start_of_table_line and line_number < end_of_table:
            read_out_genes.append(line[0])
            # print line
            # print line_number
            data_matrix.append([helper_function(value) for value in line[1:]])
            ## one line is all experimental conditions for a single gene

data_matrix = np.array(data_matrix)

# check that we got everything properly
logging.info('deleted genes:\t%s', len(deleted_genes))
logging.info('genes in read-out:\t%s', len(read_out_genes))
logging.info('data matrix shape:\t%s', data_matrix.shape)

# Show overall factor distribution
data = data_matrix.flatten()
fltr = np.logical_not(np.isnan(data))
density = gaussian_kde(data[fltr].flatten())
xs = np.linspace(data[fltr].min(), data[fltr].max(), 200)
plt.title('distribution of gene expression alteration in response to TF deletion')
ax = plt.plot(xs, -np.sqrt(-np.log(density(xs))), 'k')
plt.plot(xs, -np.sqrt(-np.log(norm.pdf(xs))), 'r')
plt.xlabel('z-score deviation')
plt.ylabel('distribution density (log)')
plt.show()


data = np.power(data_matrix, 2)
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
data2 = np.sqrt(np.nanmean(np.power(data_matrix, 2), axis=1))
## per experimental condition
fltr = np.logical_not(np.isnan(data2))
density = gaussian_kde(data2[fltr].flatten())
xs = np.linspace(data2[fltr].min(), data2[fltr].max(), 200)
plt.title('distribution of euclidean distances of gene expression deviation\nper-condition')
ax = plt.semilogy(xs, density(xs), 'k')
plt.xlabel('mean L2 z-score deviation')
plt.ylabel('distribution density (log)')
plt.show()

data3 = np.sqrt(np.nanmean(np.power(data_matrix, 2), axis=0))
## per gene
fltr = np.logical_not(np.isnan(data3))
density = gaussian_kde(data3[fltr].flatten())
xs = np.linspace(data3[fltr].min(), data3[fltr].max(), 200)
plt.title('distribution of euclidean distances of gene expression deviation\nper-gene')
ax = plt.semilogy(xs, density(xs), 'k')
plt.xlabel('mean L2 z-score deviation')
plt.ylabel('distribution density (log)')
plt.show()

# pickle.dump(data_matrix[np.logical_or(data_matrix < -2, data_matrix > 2)], open('activations.dmp', 'wb'))
pickle.dump(data_matrix, open('activations.dmp', 'wb'))
