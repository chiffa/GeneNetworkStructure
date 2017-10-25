from csv import reader as csv_reader
import numpy as np
from scipy.stats import gaussian_kde, norm, lognorm
from matplotlib import pyplot as plt

import pickle

data_rooster = []
abundance_list = []
counts_list = []

with open('nature09529-s2.csv') as source:
    reader = csv_reader(source)
    for i, line in enumerate(reader):

        if i == 20:
            print line[11:-2]
            for j, string in enumerate(line):
                if 'dNSAF AVG' in string:
                    data_rooster.append((j, j+1))

        if i > 20 and i < 2777:
            print i
            loc_abundance = []
            loc_counts = []
            for average, counts in data_rooster:
                loc_abundance.append(line[average])
                loc_counts.append(line[counts])
            abundance_list.append(loc_abundance)
            counts_list.append(loc_counts)

abundance_arr = np.array(abundance_list).astype(np.float)
counts_arr = np.array(counts_list).astype(np.int)

sel_mask_1 = counts_arr[:, 0] > 1
counts_arr = counts_arr[sel_mask_1, :]
abundance_arr = abundance_arr[sel_mask_1, :]

norm_abundance_arr = abundance_arr / np.expand_dims(abundance_arr[:, 0], 1)

print norm_abundance_arr
print counts_arr

data = norm_abundance_arr[:, 1:]
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

pickle.dump(norm_abundance_arr, open('norm_aneup.dmp', 'wb'))
