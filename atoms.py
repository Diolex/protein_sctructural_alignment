import numpy as np
import sys


def parsePdb(path):
	xyz = []
	try:
		# with open(path) as f:
			# for line in f:
				# l = line.split()
				# if l[0] == 'ATOM':
					# xyz.append([np.float32(i) for i in l[6:9]])

		#open the file
		pdbFile = open(path, 'rt')
		
		#get the name of protein
		line = pdbFile.readline()
		l = line.split()
		name = l[-1]
		# print (name)
		line = pdbFile.readline()
		
		while line:
			l = line.split()
			# skip non-CA atoms
			if not ( (l[0] == 'ATOM') and l[2] == "CA" and (line[16] == " " or line[16] == "A" ) ):
				line = pdbFile.readline();
				continue;

			xyz.append([np.float32(i) for i in l[6:9]])
			line = pdbFile.readline()
			
	except IOError:
		print("Could not open file '%s'" % self.path)
	return xyz, name

class Atoms:
	def __init__(self, xyz, name = "default"):
		self.name = name
		self.xyz = np.matrix(xyz)
		self.count = len(self.xyz)

		# if len(self.xyz) > 1:
		# 	print("x\t\ty\t\tz")
		# 	print(self.xyz)
		# 	print("Total atoms = %d" % self.count)
		# 	print (name)
		
	def centroid(self):
		return self.xyz.mean(0)

	def centerOrigin(self):
		self.xyz = np.subtract(self.xyz, self.centroid())

	def translate(self, t):
		self.xyz = np.add(self.xyz, t)
		
if __name__ == "__main__":
	if len(sys.argv) > 1:
		path = sys.argv[1]
		xyz, name = parsePdb(path)
		atoms = Atoms(xyz, name)
		print atoms.name