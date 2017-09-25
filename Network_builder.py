import pickle
import random
from itertools import product, combinations
from scipy.stats import gaussian_kde, norm
from matplotlib import pyplot as plt

import numpy as np

essentiality = 0.3  # epsilon - if pathway throughput falls below this value we consider pathway to fail
aneup_essentiality = 0.4  # espilon for aneuploid pathway
min_length = 1.99  # minimal length of the pathway
total_genes = 6500  # total genes in yeast
independent_pathways = 55  # estimation of independent pathways in yeast

length_width = pickle.load(open('w_l_accumulator.dmp', 'rb'))
activations = pickle.load(open('activations.dmp', 'rb'))
aneup_perturbation = pickle.load(open('norm_aneup.dmp', 'rb'))

aneup_perturbation = aneup_perturbation.flatten()

length_width = np.array(length_width)
_length = np.array(length_width)[:, 0]

_width = np.array(length_width)[:, 1]
_width = _width[_length > min_length]
_length = _length[_length > min_length]

length_width = np.column_stack((_length, _width))

activations = np.array(activations).flatten()
activations = activations[np.logical_not(np.isnan(activations))]

# checking if gaussian hump plays a role:
# activations[np.logical_and(activations < 1.0, activations > -1.0)] = 0  # Nope, percolation....
# activations = activations[np.logical_or(activations > 2.0, activations < -2.0)]

def simulation_run():

    essential_genes = 0
    total_genes = 0

    synthetic_lethals = 0
    total_interactions = 0

    cond_essentail_genes = 0

    for _ in range(0, 100):
        _l, _w = random.choice(length_width.tolist())
        _l = int(round(_l))
        _w = int(round(_w))
        matrix_chain = []
        abs_m_chain = []

        for _ in range(0, _l):
            mat = np.random.choice(activations, size=(_w, _w))
            matrix_chain.append(mat)  # meh - sparsity is bad
            abs_m_chain.append(np.abs(mat))

        signal = np.ones((_w, 1)) / float(_w)
        for mat in abs_m_chain:
            signal = np.dot(mat, signal)
        unperturbed = np.sum(signal)

        essentials = []
        non_essentials = []

        for i, j in product(range(0, _w), range(0, _l)):
            pad = np.ones((_l, _w))
            pad[j, i] = 0

            signal = np.ones((_w, 1)) / float(_w)
            for k, mat in enumerate(abs_m_chain):
                signal = np.dot(mat, signal)
                signal = signal * np.expand_dims(pad[k, :], 1)
            perturbed = np.sum(signal)
            if perturbed / unperturbed < essentiality:
                essential_genes += 1
                essentials.append([j, i])
            else:
                non_essentials.append([j, i])

        total_genes += _l * _w

        for j, i in essentials:
            pad = np.ones((_l, _w))
            pad[j, i] = 0

            # signal = np.expand_dims(np.array(lognorm.rvs(0.25, size=_w))/float(_w), 1)
            signal = np.ones((_w, 1)) / float(_w)
            aneuploid = np.ones((_w, 1)) / float(_w)
            no_aneup = np.ones((_w, 1)) / float(_w)
            for k, mat in enumerate(abs_m_chain):
                signal = np.dot(mat, signal)
                aneuploid = np.dot(mat, aneuploid)
                no_aneup = np.dot(mat, no_aneup)
                aneuploid_correctives = np.random.choice(aneup_perturbation, size=(_w, 1))
                aneuploid = aneuploid * aneuploid_correctives
                signal = signal * aneuploid_correctives
                signal = signal * np.expand_dims(pad[k, :], 1)
                no_aneup = no_aneup * np.expand_dims(pad[k, :], 1)
            perturbed = np.sum(signal)
            aneuploid = np.sum(aneuploid)
            no_aneup = np.sum(no_aneup)
            if perturbed / unperturbed < aneup_essentiality:  # we assume a bit of additional edge needed to overcome aneup fitness penalty
                pass
            else:
                # print '>>>>'
                # print unperturbed
                # print no_aneup
                # print aneuploid
                # print perturbed
                # print '<<<<'
                cond_essentail_genes += 1

        for gene_1, gene_2 in combinations(non_essentials, 2):
            pad = np.ones((_l, _w))
            pad[gene_1[0], gene_1[1]] = 0
            pad[gene_2[0], gene_2[1]] = 0

            signal = np.ones((_w, 1)) / float(_w)
            for k, mat in enumerate(abs_m_chain):
                signal = np.dot(mat, signal)
                signal = signal * np.expand_dims(pad[k, :], 1)
            perturbed = np.sum(signal)
            if perturbed / unperturbed < essentiality:
                synthetic_lethals += 1

            total_interactions += 1

        total_interactions += len(non_essentials) * (5500 / independent_pathways)  # basically the proteins they are not in the same pathway

    print 'current score:', float(essential_genes) / float(total_genes)
    print 'synth lethals:', float(synthetic_lethals) / float(total_interactions)
    print 'essentials no more after aneuploidy', float(cond_essentail_genes)/float(essential_genes)
    return float(essential_genes) / float(total_genes),\
           float(synthetic_lethals) / float(total_interactions),\
           float(cond_essentail_genes)/float(essential_genes)


def improved_plot(data, stats_of_interest, x_axis):
    mean = np.nanmean(data)
    std = np.nanstd(data)
    print '%s; mean: %.4f, std: %.4f' % (stats_of_interest, mean, std)
    fltr = np.logical_not(np.isnan(data))
    density = gaussian_kde(data[fltr].flatten())
    xs = np.linspace(data[fltr].min(), data[fltr].max(), 200)
    plt.title('Distribution of %s' % stats_of_interest)
    plt.plot(xs, density(xs), 'k', label='Distribution by simulation')
    plt.plot(xs, norm.pdf(xs, mean, std), 'r', label='Normal fitted to the simulation')
    maxh = np.max(norm.pdf(xs, mean, std))
    plt.vlines([mean - 2*std, mean, mean + 2*std], 0, maxh)
    plt.xlabel(x_axis)
    plt.ylabel('distribution density')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    essentials_accumulator = []
    lethal_accumulator = []
    conditional_essentials = []

    for _ in range(0, 25):
        essentials_fraction, lethal_int_fraction, cond_ess_fraction = simulation_run()
        essentials_accumulator.append(essentials_fraction)
        lethal_accumulator.append(lethal_int_fraction)
        conditional_essentials.append(cond_ess_fraction)

    improved_plot(np.array(essentials_accumulator)*100,
                  'essential genes prevalence',
                  'percentage of all genes')
    improved_plot(np.array(lethal_accumulator)*100,
                  'lethal interactions prevalence',
                  'percentage of all non-essential gene pairs')
    improved_plot(np.array(conditional_essentials)*100,
                  'evolvable essential genes',
                  'percentage of all essential genes')


