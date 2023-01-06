import treeswift as ts
from math import exp
import numpy as np
import multiprocessing as mp
from sys import platform as _platform
from sklearn.linear_model import LinearRegression

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

def set_refs(reference, t_dist):
    for i, ref in enumerate(reference.refs):
        yield (i, ref, t_dist[ref], reference.refs,)


def get_weights_all(reference, tree, options):
    print("Calculating weights using LR ...")

    # get the dist matrix from tree
    t_dist = tree.distance_matrix(True)
    refs = reference.refs
    L = len(list(refs.items())[0][1])
    n = len(refs)
    X = np.zeros((int(n*(n-1)/2), L))
    y = np.zeros((int(n*(n-1)/2), ))
    cnt = 0

    for i, ref1 in enumerate(refs):
        a2 = refs[ref1]
        for j, ref2 in enumerate(refs):
            if i < j:
                b2 = refs[ref2]
                nondash = np.logical_and(a2 != b'-', b2 != b'-')
                X[cnt] = np.logical_and(a2 != b2, nondash)
                y[cnt] = jc_inverse(t_dist[ref1][ref2])
                cnt += 1
        print(i, ref1, "done")

    # dataset prepared
    lr = LinearRegression(fit_intercept=False)

    lr.fit(X, y)
    w = lr.coef_
    
    if options.weight_output_fp:
        with open(options.weight_output_fp, "w+") as fp:
            fp.write(" ".join((str(x) for x in w)))

    return w


def get_weights(reference, tree, options):
    print("Calculating weights ...")

    # get the dist matrix from tree
    t_dist = tree.distance_matrix(True)
    refs = reference.refs
    L = len(list(refs.items())[0][1])

    w = np.zeros(L)
    denom = np.ones(L)
    num = np.zeros(L)

    if _platform == "win32" or _platform == "win64" or _platform == "msys":
        # if windows, multithreading is not supported until either
        # processes can be forked in windows or apples works with spawn.
        results_combined = list(map(lambda a: calc(a[0], *a[1:]), set_refs(reference, t_dist)))   # a.k.a starmap
    else:
        pool = mp.Pool(options.num_thread)
        combined = pool.starmap(calc, set_refs(reference, t_dist))
    
    for c in combined:
        num += c[0]
        denom += c[1]

    w = num / denom
    w /= np.sum(w)
    
    print("Weights calculated.")

    if options.weight_output_fp:
        with open(options.weight_output_fp, "w+") as fp:
            fp.write(" ".join((str(x) for x in w)))

    optimization_results(reference.refs, t_dist, w)

    return w 


def optimization_results(refs, t_dist, w):
    print("Before and after...")
    total_before = 0; total_after = 0
    w_b = np.ones(len(w)) / len(w)
    for i, ref in enumerate(refs):
        a2 = refs[ref]
        for j, ref2 in enumerate(refs):
            if i > j:
                b2 = refs[ref2]
                # total_before += 
                nondash = np.logical_and(a2 != b'-', b2 != b'-')
                hh = np.logical_and(a2 != b2, nondash)
                tt = jc_inverse(t_dist[ref][ref2])

                total_before += np.sum(np.square(w_b*hh - tt))
                total_after += np.sum(np.square(w*hh - tt))   

    print(total_before, total_after)
    print("Decrease = ", (total_before-total_after)*100.0/total_before, "percent")
    exit(0)


# def get_weights_old(reference, tree, options):
#     print("Calculating weights ...")

#     # get the dist matrix from tree
#     t_dist = tree.distance_matrix(True)
#     refs = reference.refs
#     L = len(list(refs.items())[0][1])

#     w = np.zeros(L)
#     denom = np.ones(L)
#     num = np.zeros(L)

#     for i, ref in enumerate(refs):
#         a2 = refs[ref]
#         for j, ref2 in enumerate(refs):
#             if j > i:
#                 b2 = refs[ref2]
#                 nondash = np.logical_and(a2 != b'-', b2 != b'-')
#                 hh = np.logical_and(a2 != b2, nondash)
#                 tt = jc_inverse(t_dist[ref][ref2])
#                 num += tt * hh
#                 denom += hh

#     w = num / denom  
#     w /= np.sum(w)
    
#     print("Weights calculated.")

#     if options.weight_output_fp:
#         with open(options.weight_output_fp, "w+") as fp:
#             fp.write(" ".join((str(x) for x in w)))

#     return w 


def calc(i, ref, t_dist, refs):
    a2 = refs[ref]
    L = len(a2)
    denom = np.zeros(L)
    num = np.zeros(L)
    for j, ref2 in enumerate(refs):
        if j > i:
            b2 = refs[ref2]
            nondash = np.logical_and(a2 != b'-', b2 != b'-')
            hh = np.logical_and(a2 != b2, nondash)
            tt = jc_inverse(t_dist[ref2])
            num += tt * hh
            denom += hh

    return num, denom
    