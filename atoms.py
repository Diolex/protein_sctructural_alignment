import numpy as np
import sys


def parsePdb(path):
    xyz = []
    try:
        with open(path) as f:
            for line in f:
                l = line.split()
                if l[0] == 'ATOM':
                    xyz.append([np.float32(i) for i in l[6:9]])
    except IOError:
        print("Could not open file '%s'" % self.path)
    return xyz

class Atoms:
    def __init__(self, atoms):
        self.xyz = np.matrix(atoms)
        self.count = len(self.xyz)
        #if len(self.xyz) > 1:
        #    print("x\t\ty\t\tz")
        #    print(self.xyz)
        #    print("Total atoms = %d" % self.count)

    def centroid(self):
        return self.xyz.mean(0)

    def centerOrigin(self):
        self.xyz = np.subtract(self.xyz, self.centroid())

    def translate(self, t):
        self.xyz = np.add(self.xyz, t)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        path = sys.argv[1]
        atoms = Atoms(parsePdb(path))
