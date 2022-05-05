import numpy as np
import math
import sys
from matplotlib import pyplot as plt
from matplotlib import colors

def main(mapFile, outFile):
    mapMatrix = []

    dim = 0
    nproc = 0
    length = 0
    with open(mapFile, 'r') as f:
        line = f.readline()
        dim, dim = line.split()
        dim = int(dim)

        for i in range(dim):
            line = f.readline()
            mapMatrix.append(line.split())

        mapMatrix = [[int(x) for x in row] for row in mapMatrix]

    with open(outFile, 'r') as f:
        line = f.readline()
        dim, dim, nproc = line.split()
        nproc = int(nproc)
        dim = int(dim)

        for t in range(nproc):
            for x in range(dim):
                line = f.readline()
                for y, val in enumerate(line.split()): 
                    if (int(val)):
                        mapMatrix[x][y] = int(2)

        while True:
            line = f.readline()
            if not line:
                break
            x, y = line.split()
            if (x == 'LENGTH') :
                length = int(y)
            else:
                x = int(x)
                y = int(y)
                mapMatrix[x][y] = int(3)

    # mapMatrix[0][0] = int(0)
    # print(mapMatrix)

    # Color Grid with map
    cmap = colors.ListedColormap(['White','Black','Red', 'Green'])
    plt.figure(figsize=(dim,dim))
    plt.pcolor(mapMatrix[::-1],cmap=cmap,edgecolors='k', linewidths=3)
    plt.xlabel('Length: {}'.format(length))
    plt.show()

if __name__ == "__main__":
    if (len(sys.argv) != 3):
        print('Usage: python3 visualizerMPI.py <inputFile (map)> <outputFilename (results)>')
        exit
    mapFile = str(sys.argv[1])
    outFile = str(sys.argv[2])
    main(mapFile, outFile)