import numpy as np
from csv import reader as csv_reader
from matplotlib import pyplot as plt


input_file = 'essential_in_conditions.tsv'

data_array = []

with open(input_file, 'r') as input:
    reader = csv_reader(input)
    for line in reader:
        data_array.append(line)

data_array = np.array(data_array).astype(np.float)

plt.errorbar(data_array[:, 0], data_array[:, 1], yerr=data_array[:, 2], color='k', linewidth=2., elinewidth=1.)
plt.title('Environment-specific essential genes')
plt.xlabel('Environments')
plt.ylabel('Percentage of genes essential in at least one environment')
plt.xlim(0, 1010)
plt.show()
