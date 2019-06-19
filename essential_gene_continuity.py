import numpy as np
from csv import reader as csv_reader
from matplotlib import pyplot as plt
import logging
import sys
from scipy.stats import gaussian_kde, norm, chi2
import pickle

input_file = 'essentiality_threshold-ranks.tsv'

headers = []
data_matrix = []

with open(input_file, 'r') as in_file:
    reader = csv_reader(in_file, delimiter='\t')
    headers += reader.next()
    for line_number, line in enumerate(reader):
        # print line_number, line
        data_matrix.append(float(line[-1]))

data_matrix = np.array(data_matrix)

data = data_matrix.flatten()
fltr = np.logical_not(np.isnan(data))
density = gaussian_kde(data[fltr].flatten())
xs = np.linspace(data[fltr].min(), data[fltr].max(), 200)
plt.title('distribution of RANKS scores for gene essentiality_threshold')
ax = plt.plot(xs, -np.sqrt(-np.log(density(xs))), 'k')
plt.vlines([-2], 0, -3, 'r')
# plt.plot(xs, -np.sqrt(-np.log(norm.pdf(xs))), 'r')
plt.xlabel('RANKS Score')
plt.ylabel('distribution density (log)')
plt.show()