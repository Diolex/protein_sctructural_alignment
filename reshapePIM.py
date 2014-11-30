import csv

c = []
with open('results.pim', 'r') as pim:
	names = []
	for line in pim:
		l = line.split()
		if len(l) > 15:
			names.append(l[1])
			c.append(l[2:len(l)])
			
with open('wRMSD.csv', 'rb') as input:
	with open('compare.csv', 'wb') as csvfile:
		fieldnames = ['Protein A', 'Protein B','wRMSD', 'PI']
		reader = csv.reader(input)
		writer = csv.DictWriter(csvfile, fieldnames=fieldnames)	
		for row in reader:
			name_A = row[0]
			name_B = row[1]
			wRMSD = row[2]
			i = names.index(name_A)
			j = names.index(name_B)
			writer.writerow({'Protein A': name_A, 'Protein B': name_B, 'wRMSD':wRMSD,  'PI': c[i][j]})

