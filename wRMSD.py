import numpy as np
import sys
from atoms import Atoms
from atoms import parsePdb
from horns import horn

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Enter 2 .pdb files to compare")
        sys.exit()

    # get the 2 proteins loaded
    pro = [Atoms(parsePdb(sys.argv[1])), Atoms(parsePdb(sys.argv[2]))]
    mmin = min(pro[0].count, pro[1].count)
    ### STEP 1 ###
    for p in pro: # resize
        p.count = mmin
        p.xyz = np.resize(p.xyz, (mmin, 3)) # cut off extra atoms (shh)
        p.centerOrigin() # move centroid to Origin

    ### STEP 2 ###
    # average Structure
    avgS = Atoms(pro[0].xyz*0) # initialize all 0s
    for p in pro:
        avgS.xyz += p.xyz
    avgS.xyz = avgS.xyz/len(pro) # average

    w = 1 # weights not used yet
    # calculate the Standard Deviation
    sd = np.float64(0)
    for p in pro:
        for atom in range(len(p.xyz)):
            n = np.linalg.norm(p.xyz[atom] - avgS.xyz[atom])
            sd += w * n * n
    
    # e is epsilon aka reaaaally small number
    eps = 1.0*np.power(10.0,-5.0)
    ### STEP 3 ###
    while True:
        for p in pro:
            # Ri = from Horn's method
            # Si = Ri * Si
            w = np.array([1.0]*p.count)
            r,t = horn(avgS.xyz, p.xyz, w)
            p.xyz = p.xyz*r + t

        ### STEP 4 ###
        # average Structure
        avgS.xyz = avgS.xyz*0 # initialize all 0s
        for p in pro:
            avgS.xyz += p.xyz
        avgS.xyz = avgS.xyz/len(pro) # average

        # calculate the Standard Deviation
        sdNew = np.float64(0)
        for p in pro:
            for atom in range(len(p.xyz)):
                n = np.linalg.norm(p.xyz[atom] - avgS.xyz[atom])
                sdNew += w * n * n
    
        ### STEP 5 ###
        if (sd - sdNew) < eps:
            break # algorithm terminates
        else:
            sd = sdNew
