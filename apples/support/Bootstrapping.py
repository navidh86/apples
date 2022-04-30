import numpy as np
import os
from os.path import dirname
import subprocess as sp
import json
import tempfile
import shutil

class Bootstrapping:
    boot = [] # the 2d matrix
    boot2 = [] # 3d version of the same matrix
    sample_count = 0
    sequence_length = 0
    np.random.seed(56)

    @classmethod
    def get_boot_matrix(cls, sample_count, sequence_length):
        # create bootstrap matrix
        mat = np.random.choice(sequence_length, (sample_count, sequence_length), replace=True)

        # count of each item per row
        Bootstrapping.boot = np.zeros((sample_count+1, sequence_length))
        Bootstrapping.boot2 = np.zeros((sample_count+1, 1, sequence_length))
        
        # first row is all ones
        Bootstrapping.boot[0] = np.ones(sequence_length)
        Bootstrapping.boot2[0] = np.ones((1, sequence_length))

        for i in range(1, sample_count+1):
            Bootstrapping.boot[i] = np.bincount(mat[i-1], minlength=sequence_length)
            Bootstrapping.boot2[i] = Bootstrapping.boot[i].reshape((1, sequence_length))

        return Bootstrapping.boot

    @classmethod
    def perform_slow_bootstrapping(cls, tree_fp, refs, queries, sample_count, sequence_length, 
                                old_results, execpath, options):
        results = []
        for result in old_results:
            results.append({0: result})

        # refs = {}
        # for ref in old_refs:
        #     refs[ref] = old_refs[ref].tostring().decode('utf-8')
        
        # queries = {}
        # for query in old_queries:
        #     queries[query] = old_queries[query].tostring().decode('utf-8')

        for i in range(sample_count):
            print("Bootstrap sample: " + str(i+1))

            boot = np.random.choice(sequence_length, sequence_length, replace=True)
            boot_refs = {}
            boot_queries = {}

            for ref in refs:
                boot_refs[ref] = refs[ref][boot].tostring().decode('utf-8')
            for query in queries:
                boot_queries[query] = queries[query][boot].tostring().decode('utf-8')
            
            ref_fp = tempfile.NamedTemporaryFile(mode='w+t', delete=False)
            for ref in boot_refs:
                ref_fp.write(">" + ref + "\n")
                ref_fp.write(''.join(boot_refs[ref]) + "\n")
            ref_fp.close()

            query_fp = tempfile.NamedTemporaryFile(mode='w+t', delete=False)
            for query in boot_queries:
                query_fp.write(">" + query + "\n")
                query_fp.write(''.join(boot_queries[query]) + "\n")
            query_fp.close()

            out_fp = tempfile.NamedTemporaryFile(delete=False)
            out_fp.close()

            # Run APPLES on the sampled query and ref alignments
            command = ['python3', execpath, '-t', tree_fp, '-q', query_fp.name, '-s', ref_fp.name, '-o', out_fp.name, '-T', str(options.num_thread),
                    '-b', str(options.base_observation_threshold), '-f', str(options.filt_threshold), '-V', str(options.minimum_alignment_overlap)]
            if options.mask_lowconfidence:
                command.append('-X')
            if options.exclude_intplace:
                command.append('--exclude')
            if options.protein_seqs:
                command.append('-p')
            
            sp.run(command, stdout=sp.DEVNULL, stderr=sp.STDOUT)

            fp = open(out_fp.name)
            jp = fp.read()
            fp.close()
            jp = json.loads(jp)

            for idx, placement in enumerate(jp["placements"]):
                ar = [placement]
                results[idx][i+1] = {"placements": ar}

            os.remove(ref_fp.name)
            os.remove(query_fp.name)
            os.remove(out_fp.name)

        return results
