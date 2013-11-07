import csv
f = open("peptide-catalysis/recommendation/Recommendation_Sep_27/recommendation.txt")
charge = []
for peptide in f:
	onecharge = 0
	for AA in peptide:
		if (AA == 'K' or AA == 'R'):
			onecharge += 1
		elif (AA == 'E' or AA == 'D'):
			onecharge -= 1
	charge.append(onecharge)
f.close()

out = []
for element in charge:
	out.append([element, element])
with open('output.csv','wb') as ff:
	writer = csv.writer(ff)
	writer.writerows(out)
