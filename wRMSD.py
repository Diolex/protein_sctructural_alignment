import numpy as np
import sys
from atoms import Atoms
from atoms import parsePdb
from horns import horn
import math

def filter (atom_1, atom_2, ali):
	seq_1 = ""
	seq_2 = ""

	xyz_1 = []
	xyz_2 = []

	# extract the proteins
	for line in ali:
		if (line.startswith(atom_1.name) and line[5] == 'A'):
			l = line.split()
			seq_1 += l[1]
		if (line.startswith(atom_2.name) and line[5] == 'A'):
			l = line.split()
			seq_2 += l[1]			

	count_1 = 0
	count_2 = 0
	mask_1 = []
	mask_2 = []
	
	# find aligned atoms
	for i in range(len(seq_1)):
		if (seq_1[i] != '-'):
			if (seq_2[i] != '-'):
				mask_1.append(count_1)
				mask_2.append(count_2)
				count_2 += 1
			count_1 += 1
		else:
			if (seq_2[i] != '-'):
				count_2 += 1
	# print (mask_1)
	# print (mask_2)
	
	# generate new xyz
	xyz_1 = atom_1.xyz[mask_1]
	xyz_2 = atom_2.xyz[mask_2]
	
	return Atoms(xyz_1), Atoms(xyz_2)

if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Enter 2 .pdb files to compare")
		sys.exit()


	# load alignment
	try: 
		ali = open(sys.argv[3], 'rt')
	except IOError:
		print("Could not open file '%s'" % sys.argv[3])
		
	# get the 2 proteins loaded
	xyz, name = parsePdb(sys.argv[1])
	atom_1 = Atoms(xyz, name)
	xyz, name = parsePdb(sys.argv[2])	
	atom_2 = Atoms(xyz, name)
	
	# filter non-aligned atoms
	atom_1, atom_2 = filter (atom_1, atom_2, ali)
	
	pro = [atom_1, atom_2]
	mmin = min(pro[0].count, pro[1].count)
	### STEP 1 ###
	for p in pro: # resize
		p.count = mmin
		p.xyz = np.resize(p.xyz, (mmin, 3)) # cut off extra atoms (shh)
		p.centerOrigin() # move centroid to Origin

	### STEP 2 ###
	# average Structure
	avgS = Atoms(np.zeros(pro[0].xyz.shape)) # initialize all 0s
	for p in pro:
		avgS.xyz += p.xyz
	avgS.xyz = avgS.xyz/len(pro) # average

	w = 1 # weights not used yet
	weight = np.array([1.0]*p.count)
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
			r,t = horn(p.xyz, avgS.xyz, weight)
			#p.xyz = p.xyz*r + t
			temp = np.outer(np.ones(avgS.xyz.shape[0]), t)
			p.xyz = np.dot(p.xyz,r.T) + temp
			
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
		#print (sd)
	wrmsd = math.sqrt(sd / len(pro))
	print 'wRMSD is',wrmsd