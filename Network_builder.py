import pickle
import random
from itertools import product, combinations

import numpy as np


length_width = pickle.load(open('w_l_accumulator.dmp', 'rb'))
activations = pickle.load(open('activations.dmp', 'rb'))



length_width = np.array(length_width)
_length = np.array(length_width)[:, 0]

_width = np.array(length_width)[:, 1]
_width = _width[_length > 1.99]
_length = _length[_length > 1.99]

length_width = np.column_stack((_length, _width))

activations = np.array(activations).flatten()
activations = activations[np.logical_not(np.isnan(activations))]

# checking if gaussian hump plays a role:
activations[np.logical_and(activations < 1.0, activations > -1.0)] = 0

essential_genes = 0
total_genes = 0

synthetic_lethals = 0
total_interactions = 0

for _ in range(0, 100):
    _l, _w = random.choice(length_width.tolist())
    _l = int(round(_l))
    _w = int(round(_w))
    # print _l, _w
    matrix_chain = []
    abs_m_chain = []

    for _ in range(0, _l):
        mat = np.random.choice(activations, size=(_w, _w))
        matrix_chain.append(mat)  # meh - sparsity is bad
        abs_m_chain.append(np.abs(mat))

    signal = np.ones((_w, 1)) / float(_w)
    for mat in abs_m_chain:
        signal = np.dot(mat, signal)
        # print signal
    unperturbed = np.sum(signal)
    # print unperturbed

    non_essentials = []
    for i, j in product(range(0, _w), range(0, _l)):
        pad = np.ones((_l, _w))
        pad[j, i] = 0

        signal = np.ones((_w, 1)) / float(_w)
        for k, mat in enumerate(abs_m_chain):
            signal = np.dot(mat, signal)
            # print signal
            # print signal.shape, pad[k, :].shape
            signal = signal * np.expand_dims(pad[k, :], 1)
            # print signal.shape
            # print signal
        perturbed = np.sum(signal)
        # print "perturbed", perturbed
        # print "ratio", perturbed/unperturbed
        if perturbed/unperturbed < 0.4:
            # print "candidate gene detected"
            essential_genes += 1
        else:
            non_essentials.append([j, i])

    total_genes += _l * _w

    for gene_1, gene_2 in combinations(non_essentials, 2):
        pad = np.ones((_l, _w))
        pad[gene_1[0], gene_1[1]] = 0
        pad[gene_2[0], gene_2[1]] = 0

        signal = np.ones((_w, 1)) / float(_w)
        for k, mat in enumerate(abs_m_chain):
            signal = np.dot(mat, signal)
            # print signal
            # print signal.shape, pad[k, :].shape
            signal = signal * np.expand_dims(pad[k, :], 1)
            # print signal.shape
            # print signal
        perturbed = np.sum(signal)
        # print "perturbed", perturbed
        # print "ratio", perturbed/unperturbed
        if perturbed/unperturbed < 0.4:
            synthetic_lethals += 1

        total_interactions += 1

    total_interactions += len(non_essentials)*(5500/50) # basically the proteins they are not in the same pathway

    print 'current score:', float(essential_genes) / float(total_genes)
    print 'synth lethals:', float(synthetic_lethals) / float(total_interactions)
    # raw_input('press any key')

