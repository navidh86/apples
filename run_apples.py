#!/usr/bin/env python3
import re
from apples.PoolQueryWorker import PoolQueryWorker
from apples.fasta2dic import fasta2dic
from apples.Reference import ReducedReference
from apples.OptionsRun import options_config
import multiprocessing as mp
from apples.jutil import join_jplace, join_jplace_support, join_jplace_support_all
import sys
import json
from sys import platform as _platform

from apples.prepareTree import prepareTree
import time
import logging
import pickle

from apples.support.Bootstrapping import Bootstrapping
from apples.weight import get_weights, get_weights_slow
import numpy as np
import os.path


if __name__ == "__main__":
    mp.set_start_method('fork')
    startb = time.time()
    options, args = options_config()
    logging.info("[%s] Options are parsed." % time.strftime("%H:%M:%S"))

    if options.database_fp:
        # unpickle tree from database
        start = time.time()
        fdtb = open(options.database_fp, "rb")
        up = pickle.Unpickler(fdtb)
        first_read_tree = up.load()
        name_to_node_map = up.load()
        extended_newick_string = up.load()
        logging.info(
            "[%s] Tree is loaded from APPLES database in %.3f seconds." % (
                time.strftime("%H:%M:%S"), (time.time() - start)))

    if options.tree_fp:
        orig_tree_fp = options.tree_fp
        first_read_tree, name_to_node_map, extended_newick_string = prepareTree(options)

    if options.dist_fp:
        reference = None

        def read_dismat(f):
            tags = list(re.split("\s+", f.readline().rstrip()))[1:]
            for line in f.readlines():
                dists = list(re.split("\s+", line.strip()))
                query_name = dists[0]
                obs_dist = dict(zip(tags, map(float, dists[1:])))
                yield (query_name, None, obs_dist)

        start = time.time()
        f = open(options.dist_fp)
        queries = read_dismat(f)
        logging.info(
            "[%s] Query sequences are prepared in %.3f seconds." % (time.strftime("%H:%M:%S"), (time.time() - start)))
    else:
        if options.ref_fp:
            start = time.time()
            weighted_flag = True
            load_weights = True
            reference = ReducedReference(options.ref_fp, options.protein_seqs, options.tree_fp,
                                         options.filt_threshold, options.num_thread, weighted_flag)
            logging.info(
                "[%s] Reduced reference is computed in %.3f seconds." % (
                    time.strftime("%H:%M:%S"), (time.time() - start)))
        else:  # options.database_fp
            start = time.time()
            reference = up.load()
            fdtb.close()
            logging.info(
                "[%s] Reduced reference is loaded from APPLES database in %.3f seconds." % (
                    time.strftime("%H:%M:%S"), (time.time() - start)))

        reference.set_baseobs(options.base_observation_threshold)
        
        if weighted_flag:
            if not load_weights:
                weights = get_weights_slow(reference, first_read_tree, options)
            else:
                # load from file
                filename = "weights.txt"
                weights = open(filename, "r").read().split(" ")
                weights = [float(w) for w in weights]
        else:
            weights = None
        reference.set_weights(weights)

        start = time.time()
        if options.query_fp:
            query_dict = fasta2dic(options.query_fp, options.protein_seqs, options.mask_lowconfidence)
        else:  # must be extended reference
            extended_dict = fasta2dic(options.extended_ref_fp, options.protein_seqs, options.mask_lowconfidence)
            query_dict = {key: value for key, value in extended_dict.items() if key not in reference.refs}

        def set_queries(query_dict):
            for query_name, query_seq in query_dict.items():
                yield (query_name, query_seq, None)

        queries = set_queries(query_dict)
        logging.info(
            "[%s] Query sequences are prepared in %.3f seconds." % (time.strftime("%H:%M:%S"), (time.time() - start)))

    startq = time.time()
    queryworker = PoolQueryWorker()
    queryworker.set_class_attributes(reference, options, name_to_node_map)

    if options.find_support:
        Bootstrapping.sample_count = options.sample_count

    if not options.fast_support:
        query_function = queryworker.runquery
    else:
        boot = Bootstrapping.get_boot_matrix(options.sample_count, len(reference.representatives[0][0]))
        query_function = queryworker.runquery_support_fast

    if _platform == "win32" or _platform == "win64" or _platform == "msys":
        # if windows, multithreading is not supported until either
        # processes can be forked in windows or apples works with spawn.
        results_combined = list(map(lambda a: query_function(a[0], *a[1:]), queries))   # a.k.a starmap
    else:
        pool = mp.Pool(options.num_thread)
        results_combined = pool.starmap(query_function, queries)
    logging.info(
        "[%s] Processed all queries in %.3f seconds." % (time.strftime("%H:%M:%S"), (time.time() - startq)))

    results = {}
    valids = {}
    if options.find_support:
        if options.fast_support:
            for result in results_combined:                
                results[result[0]] = result[1]
                valids[result[0]] = result[2]
        else:
            results = []
            for result in results_combined:                
                results.append(result[0])
                valids[result[0]['placements'][0]['n'][0]] = result[1]
    else:
        results = [rc[0] for rc in results_combined]

    if not options.find_support:
        result = join_jplace(results)
    else:
        if not options.fast_support:
            results = Bootstrapping.perform_slow_bootstrapping(orig_tree_fp, reference.refs, query_dict, 
                                    options.sample_count, len(reference.representatives[0][0]), results, 
                                    os.path.abspath(__file__), options)
        result = join_jplace_support_all(results, valids, options.keep_factor, options.keep_at_most, options.prioritize_lse)                                    

    result["tree"] = extended_newick_string
    result["metadata"] = {"invocation": " ".join(sys.argv)}
    result["fields"] = ["edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length"]
    result["version"] = 3

    if options.output_fp:
        f = open(options.output_fp, "w")
    else:
        f = sys.stdout
    f.write(json.dumps(result, sort_keys=True, indent=4))
    f.write("\n")
    f.close()
    logging.warning(
        "[%s] APPLES finished in %.3f seconds." % (time.strftime("%H:%M:%S"), (time.time() - startb)))
