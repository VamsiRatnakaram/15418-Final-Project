import numpy as np
import math
import sys

def main(size, filename):
    matrix = np.zeros((size, size))
    matrix = [[(1 if np.random.rand() < 0.25 else x) for x in line] for line in matrix]

    with open(filename, 'w+') as testfile:
        testfile.write('{} {}\n'.format(size, size))
        for row in matrix:
            testfile.write(' '.join([str(int(a)) for a in row]) + '\n')

if __name__ == "__main__":
    if (len(sys.argv) != 3):
        print('Usage: python3 testGenerate <matrixDimension> <outputFilename>')
    size = int(sys.argv[1])
    filename = str(sys.argv[2])
    main(size, filename)