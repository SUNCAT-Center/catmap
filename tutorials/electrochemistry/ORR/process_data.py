f = open('100_energies.txt')
lines = f.readlines()
for line in lines:
	data = line.split(' ')
	print data[1]