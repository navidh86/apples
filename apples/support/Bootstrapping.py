import numpy as np
import os
import subprocess as sp
import json

class Bootstrapping:
    boot = [] # the 2d matrix
    boot2 = [] # 3d version of the same matrix
    sample_count = 0
    sequence_length = 0
    np.random.seed(56)

    @classmethod
    def getBootMatrix(cls, sample_count, sequence_length):
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

        Bootstrapping.sample_count = sample_count

        print("Shape = " + str(Bootstrapping.boot2.shape))

        return Bootstrapping.boot

    @classmethod
    def performSlowBootstrapping(cls, tree_fp, old_refs, old_queries, sample_count, sequence_length, old_results):
        np.random.seed(56)

        if not os.path.isdir('Boot'):
            os.mkdir('Boot')

        results = []
        for result in old_results:
            results.append({0: result})

        refs = {}
        for ref in old_refs:
            refs[ref] = old_refs[ref].tostring().decode('utf-8')
        
        queries = {}
        for query in old_queries:
            queries[query] = old_queries[query].tostring().decode('utf-8')

        for i in range(Bootstrapping.sample_count):
            boot = np.random.choice(sequence_length, sequence_length, replace=True)
            boot_refs = {}
            boot_queries = {}

            for j, pos in enumerate(boot):
                # ref
                for ref in refs:
                    if j == 0:
                        boot_refs[ref] = []
                    boot_refs[ref].append(refs[ref][pos])
                # query
                for query in queries:
                    if j == 0:
                        boot_queries[query] = []
                    boot_queries[query].append(queries[query][pos])
            
            # print(boot_queries)

            fp = open("Boot/boot_ref_"+str(i+1)+'.fa', "w")
            for ref in boot_refs:
                fp.write(">" + ref + "\n")
                fp.write(''.join(boot_refs[ref]) + "\n")
                boot_refs[ref] = []
            fp.close()

            fp = open("Boot/boot_query_"+str(i+1)+'.fa', "w")
            for query in boot_queries:
                fp.write(">" + query + "\n")
                fp.write(''.join(boot_queries[query]) + "\n")
                boot_queries[query] = []
            fp.close()

            # Run APPLES on the sample query and ref alignments
            sp.run(['python3', 'run_apples.py', '-t', tree_fp, '-q', 'Boot/boot_query_'+str(i+1)+'.fa', 
            '-s', 'Boot/boot_ref_'+str(i+1)+'.fa', '-o', 'Boot/boot_out_'+str(i+1)+'.jplace'], stdout=sp.DEVNULL, stderr=sp.STDOUT)

            fp = open("Boot/boot_out_" + str(i+1) + ".jplace")
            jp = fp.read()
            fp.close()
            jp = json.loads(jp)

            for idx, placement in enumerate(jp["placements"]):
                ar = [placement]
                results[idx][i+1] = {"placements": ar}

            os.remove('Boot/boot_query_'+str(i+1)+'.fa')
            os.remove('Boot/boot_ref_'+str(i+1)+'.fa')
            os.remove('Boot/boot_out_'+str(i+1)+'.jplace')

        os.rmdir('Boot')

        return results
