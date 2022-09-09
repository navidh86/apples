import treeswift as ts
from math import exp
import numpy as np

def jc_inverse(d):
    return 0.75 * (1 - exp(-4.0*d/3))

# def get_weights(reference, tree, options):
#     # get the dist matrix from tree
#     t_dist = tree.distance_matrix(True)
#     # print(t_dist)
#     refs = reference.refs
#     L = len(list(reference.refs.items())[0][1])
#     hh = np.zeros((len(refs), len(refs), L))
#     tt = np.zeros((len(refs), len(refs)))

#     for i, ref in enumerate(refs):
#         a2 = refs[ref]
#         for j, ref2 in enumerate(refs):
#             b2 = refs[ref2]
#             nondash = np.logical_and(a2 != b'-', b2 != b'-')
#             hh[i][j] = np.logical_and(a2 != b2, nondash)
#             if ref != ref2:
#                 tt[i][j] = jc_inverse(t_dist[ref][ref2])

#     w = np.zeros(L)
#     num = np.sum(tt[:, :, np.newaxis] * hh, axis=(0, 1))
#     denom = np.sum(hh, axis=(0, 1)) + np.ones(L)
#     # print(denom)
#     w = num / denom
    
#     w /= np.sum(w)
#     print(np.sum(w))

#     return w


def get_weights(reference, tree, options):
    print("Calculating weights ...")

    # get the dist matrix from tree
    t_dist = tree.distance_matrix(True)
    refs = reference.refs
    L = len(list(reference.refs.items())[0][1])

    w = np.zeros(L)
    denom = np.ones(L)
    num = np.zeros(L)

    for i, ref in enumerate(refs):
        a2 = refs[ref]
        for j, ref2 in enumerate(refs):
            b2 = refs[ref2]
            nondash = np.logical_and(a2 != b'-', b2 != b'-')
            hh = np.logical_and(a2 != b2, nondash)
            if ref != ref2:
                tt = jc_inverse(t_dist[ref][ref2])
            else:
                tt = 0
            num += tt * hh
            denom += hh

    w = num / denom  
    w /= np.sum(w)
    
    print("Weights calculated.")

    if options.weight_output_fp:
        with open(options.weight_output_fp, "w+") as fp:
            fp.write(" ".join((str(x) for x in w)))

    return w    