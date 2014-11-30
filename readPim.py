#check the general spread of the sequence identity of the dataset

c = [0]*10
with open('results.pim', 'r') as pim:
	for line in pim:
		l = line.split()
		if len(l) > 15:
			for n in range(3, len(l)):
				num = float(l[n])
				for m in range(0,10):
					if(num > m*10 and float(n) <= (m+1)*10):
						c[m] += 1

for m in range(0,len(c)):
	print("%.2f:%d" % (m*10.0, c[m]))
