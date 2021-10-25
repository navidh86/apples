import numpy as np

class Bootstrapping:
    boot = [] # the 2d matrix
    boot2 = [] # 3d version of the same matrix
    sample_count = 0

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