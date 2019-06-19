import pickle
import random
import numpy as np
from itertools import product, combinations
from scipy.stats import gaussian_kde, norm, chi
from matplotlib import pyplot as plt
from time import time

# declaring the variables:
essentiality_threshold = 0.3  # epsilon - if pathway throughput falls below this value we consider pathway to fail
aneup_essentiality_threshold = 0.4  # espilon for aneuploid pathway
min_length = 1.99  # minimal length of the pathway
total_non_essential_pool = 5400
independent_pathways = 17  # estimation of independent pathways in yeast
stress_conditions = 0  # how many stress conditions are we taking in account
aneuploid_passes = 20  # how many aneuploidy samples we perform to calcylate the statistics
partial_stress_essentiality = essentiality_threshold
partial_stress_essentiality = 0.6

length_width = pickle.load(open('w_l_accumulator.dmp', 'r'))
activations = pickle.load(open('activations.dmp', 'r'))
aneup_perturbation = pickle.load(open('norm_aneup.dmp', 'r'))
aneup_perturbation = aneup_perturbation.flatten()

length_width = np.array(length_width)
_length = np.array(length_width)[:, 0]

_width = np.array(length_width)[:, 1]
_width = _width[_length > min_length]
_length = _length[_length > min_length]

randomly_perturbed = False

activations = np.array(activations).flatten()
activations = activations[np.logical_not(np.isnan(activations))]

# ## uncomment if we want to have random perturbations:
# essentiality_threshold = 0.5
# aneup_essentiality_threshold = 0.6
#
# independent_pathways = 50
#
# f_length = np.random.rand(*_length.shape)
# f_length *= np.max(_length)
# f_length += 2.
#
# f_width = np.random.rand(*_width.shape)
# f_width *= np.max(_width)
# f_width += 1.
#
# _length = f_length
# _width = f_width
# #
# activations = np.random.randn(*activations.shape)
# aneup_perturbation = np.random.lognormal(0, 0.5, *aneup_perturbation.shape)
##

length_width = np.column_stack((_length, _width))

# checking if gaussian hump plays a role:
# activations[np.logical_and(activations < 1.0, activations > -1.0)] = 0  # Nope, percolation....
# activations = activations[np.logical_or(activations > 2.0, activations < -2.0)]

def simulation_run():

    essential_genes = 0
    total_genes = 0

    synthetic_lethals = 0
    total_interactions = 0

    cond_essentail_genes = 0
    genes_essential_in_cond = 0
    epistases = []

    counter = 6000

    while counter > 0:
        _l, _w = random.choice(length_width.tolist())
        _l = int(round(_l))
        _w = int(round(_w))
        counter -= _l*_w
        matrix_chain = []
        abs_m_chain = []
        single_gene_effects = {}

        for _ in range(0, _l):
            mat = np.random.choice(activations, size=(_w, _w))
            matrix_chain.append(mat)  # meh - sparsity is bad
            abs_m_chain.append(np.abs(mat))

        if randomly_perturbed:
            signal = (np.random.rand(_w, 1) + 0.5) / float(_w)
        else:
            signal = np.ones((_w, 1)) / float(_w)
        ref_signal = signal
        for mat in abs_m_chain:
            signal = np.dot(mat, signal)
        unperturbed = np.sum(signal)

        essentials = []
        non_essentials = []

        for i, j in product(range(0, _w), range(0, _l)):
            pad = np.ones((_l, _w))
            pad[j, i] = 0

            if randomly_perturbed:
                signal = (np.random.rand(_w, 1) + 0.5) / float(_w)
            else:
                signal = np.ones((_w, 1)) / float(_w)
            for k, mat in enumerate(abs_m_chain):
                signal = np.dot(mat, signal)
                signal = signal * np.expand_dims(pad[k, :], 1)
            perturbed = np.sum(signal)
            single_gene_effects[(j, i)] = perturbed / unperturbed
            if perturbed / unperturbed < essentiality_threshold:
                essential_genes += 1
                essentials.append([j, i])
            else:
                non_essentials.append([j, i])

        total_genes += _l * _w

        for j, i in non_essentials:
            pad = np.ones((_l, _w))
            pad[j, i] = 0

            flag = False
            for _ in range(0, stress_conditions):
                # signal = (np.random.beta(50, 50, (_w, 1))+0.5) / float(_w)
                # print 'stress env signal', signal/ref_signal
                signal = (np.random.rand(_w, 1) + 0.5) / float(_w)
                # signal = np.ones((_w, 1)) / float(_w)
                for k, mat in enumerate(abs_m_chain):
                    signal = np.dot(mat, signal)
                    signal = signal * np.expand_dims(pad[k, :], 1)
                perturbed = np.sum(signal)
                if perturbed / unperturbed < partial_stress_essentiality:
                    flag = True
            if flag:
                genes_essential_in_cond += 1

        for j, i in essentials:
            pad = np.ones((_l, _w))
            pad[j, i] = 0
            aneuploid_pad = []

            # signal = np.expand_dims(np.array(lognorm.rvs(0.25, size=_w))/float(_w), 1)

            conditional_essential = False
            for _ in range(0, aneuploid_passes):
                if randomly_perturbed:
                    signal = (np.random.rand(_w, 1) + 0.5) / float(_w)
                else:
                    signal = np.ones((_w, 1)) / float(_w)
                aneuploid = np.ones((_w, 1)) / float(_w)
                no_aneup = np.ones((_w, 1)) / float(_w)
                for k, mat in enumerate(abs_m_chain):
                    signal = np.dot(mat, signal)
                    aneuploid = np.dot(mat, aneuploid)
                    no_aneup = np.dot(mat, no_aneup)

                    aneuploid_correctives = np.random.choice(aneup_perturbation, size=(_w, 1))
                    aneuploid_pad.append(aneuploid_correctives)

                    aneuploid = aneuploid * aneuploid_correctives

                    signal = signal * aneuploid_correctives
                    signal = signal * np.expand_dims(pad[k, :], 1)

                    no_aneup = no_aneup * np.expand_dims(pad[k, :], 1)

                perturbed = np.sum(signal)
                aneuploid = np.sum(aneuploid)
                no_aneup = np.sum(no_aneup)
                if unperturbed / perturbed < aneup_essentiality_threshold:  # we assume a bit of additional edge needed to overcome aneup fitness penalty
                    conditional_essential = True
                #     print '----'
                #     print pad
                #     print abs_m_chain
                #     print aneuploid_pad
                #     print '>>>>'
                #     print unperturbed
                #     print no_aneup
                #     print 'f %0.3f' % (float(no_aneup)/float(unperturbed))
                #     print aneuploid
                #     print perturbed
                #     print 'f %0.3f' % (float(perturbed)/float(aneuploid))
                #     print '<<<<'
                # else:
                #     '~~~~~'
            if not conditional_essential:
                # print '!!!!!!!!'
                pass
            else:
                cond_essentail_genes += 1


        for gene_1, gene_2 in combinations(non_essentials, 2):
            pad = np.ones((_l, _w))
            pad[gene_1[0], gene_1[1]] = 0
            pad[gene_2[0], gene_2[1]] = 0

            if randomly_perturbed:
                signal = (np.random.rand(_w, 1) + 0.5) / float(_w)
            else:
                signal = np.ones((_w, 1)) / float(_w)
            for k, mat in enumerate(abs_m_chain):
                signal = np.dot(mat, signal)
                signal = signal * np.expand_dims(pad[k, :], 1)
            perturbed = np.sum(signal)
            g_eff_1 = single_gene_effects[(gene_1[0], gene_1[1])]
            g_eff_2 = single_gene_effects[(gene_2[0], gene_2[1])]
            g_eff_12 = perturbed / unperturbed
            epistases.append((g_eff_1, g_eff_2, np.log2(g_eff_12) - np.log2(g_eff_1) - np.log2(g_eff_2)))
            if perturbed / unperturbed < essentiality_threshold:
                synthetic_lethals += 1

            total_interactions += 1

        # we've calculated the genes with which they are in teh same pathway. Now we need to get a
        # sample of genes with which they proteins are not interacting
        finite_population_corrector = 1.-1./(total_non_essential_pool-len(non_essentials))
        total_interactions += len(non_essentials)*(len(non_essentials)-1)/2*(independent_pathways-1)*finite_population_corrector


        # basically the biomolecules that are not in the same pathway - that pathway representing only 1/#pathways options

        # expected number of independent genes when we poll x genes out of 5500

    epistases = np.array(epistases)

    print '\tcurrent score:', float(essential_genes) / float(total_genes)
    print '\tsynth lethals:', float(synthetic_lethals) / float(total_interactions)
    print '\tessentials no more after aneuploidy', float(cond_essentail_genes)/float(
        essential_genes)
    print '\tgenes essential in condition', float(genes_essential_in_cond)/float(total_genes -
                                                                                essential_genes)
    print '\ttotal genes covered in round', total_genes
    return float(essential_genes) / float(total_genes),\
           float(synthetic_lethals) / float(total_interactions),\
           float(cond_essentail_genes)/float(essential_genes),\
           float(genes_essential_in_cond)/float(total_genes - essential_genes),\
             epistases


def improved_plot(data, stats_of_interest, x_axis):
    mean = np.nanmean(data)
    std = np.nanstd(data)
    print '%s; mean: %.4f, std: %.4f' % (stats_of_interest, mean, std)
    fltr = np.logical_not(np.isnan(data))
    density = gaussian_kde(data[fltr].flatten())
    xs = np.linspace(data[fltr].min(), data[fltr].max(), 200)
    plt.title('Distribution of %s\n mean: %.2f; std %.2f; ' % (stats_of_interest, mean, std))
    plt.plot(xs, density(xs), 'k', label='Distribution by simulation')
    # plt.plot(xs, norm.pdf(xs, mean, std), 'r', label='Normal fitted to the simulation')
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
    ess_in_cond = []

    previous = time()
    for _i in range(0, 500):
        essentials_fraction, lethal_int_fraction, cond_ess_fraction, ess_in_cond_fraction, epistases = simulation_run()
        essentials_accumulator.append(essentials_fraction)
        lethal_accumulator.append(lethal_int_fraction)
        conditional_essentials.append(cond_ess_fraction)
        ess_in_cond.append(ess_in_cond_fraction)
        print 'round %s; %.2f s' % (_i, time() - previous)
        previous = time()

    improved_plot(np.array(essentials_accumulator)*100,
                  'essential genes prevalence',
                  'percentage of all genes')
    improved_plot(np.array(lethal_accumulator)*100,
                  'lethal interactions prevalence',
                  'percentage of all non-essential gene pairs')
    improved_plot(np.array(conditional_essentials)*100,
                  'evolvable essential genes',
                  'percentage of all essential genes')
    improved_plot(np.array(ess_in_cond)*100,
                  'condition-specific essential in at least 1 out of 250 conditions',
                  'percentage of all non-essential genes')

    plt.plot(np.sqrt(epistases[:, 0]*epistases[:, 1]), epistases[:, 2], 'ko')

    plt.show()

    plt.plot(_length, _width, 'ko')
    plt.ylabel('width of the pathway')
    plt.xlabel('length of the pathway')
    plt.show()
