import pickle
import numpy as np
from csv import reader as csv_reader
from scipy.stats import gaussian_kde, norm, lognorm
from matplotlib import pyplot as plt

average_counts_columns = []
strain_abundance_table = []
strain_counts_table = []

with open('nature09529-s2.csv') as source:
    reader = csv_reader(source)
    for line_no, line in enumerate(reader):

        # slip the header data (line_no<20)

        if line_no == 20:  # read the titles line
            print line[11:-2]
            for j, string in enumerate(line):
                if 'dNSAF AVG' in string:
                    average_counts_columns.append((j, j + 1))

        if 2777 > line_no > 20:  # read the data
            print line_no
            gene_in_strain_abundance = []
            gene_in_strain_counts = []
            for average, counts in average_counts_columns:
                gene_in_strain_abundance.append(line[average])
                gene_in_strain_counts.append(line[counts])
            strain_abundance_table.append(gene_in_strain_abundance)
            strain_counts_table.append(gene_in_strain_counts)

strain_abundance_array = np.array(strain_abundance_table).astype(np.float)
strain_counts_array = np.array(strain_counts_table).astype(np.int)

sel_mask_1 = strain_counts_array[:, 0] > 1  # we don't want to compare to 0
strain_counts_array = strain_counts_array[sel_mask_1, :]
strain_abundance_array = strain_abundance_array[sel_mask_1, :]

normalized_strain_abundance_array = strain_abundance_array / np.expand_dims(strain_abundance_array[:, 0], 1)

print normalized_strain_abundance_array
print strain_counts_array

data = normalized_strain_abundance_array[:, 1:]
fltr = np.logical_not(np.isnan(data))
density = gaussian_kde(data[fltr].flatten())
xs = np.linspace(data[fltr].min(), data[fltr].max(), 200)
plt.title('Distribution of protein dosage perturbation \n in aneuploid relative to euploid')
plt.semilogy(xs, density(xs), 'k')
plt.semilogy(xs, lognorm.pdf(xs, 0.5), 'r')
plt.axvline(1., color='k')
plt.xlabel('gene abundance relative to euploid ')
plt.ylabel('distribution density (log)')
plt.legend()
plt.show()

pickle.dump(normalized_strain_abundance_array, open('norm_aneup.dmp', 'wb'))
