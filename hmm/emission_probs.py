#!/usr/bin/env python3

sequences = []
labels = []

bind_counts = dict()
nonbinding_counts = dict()

with open('Aug-Train.dat') as f:
	lines = [line.rstrip() for line in f]
	for i in range(len(lines)):
		if i % 4 == 2:
			sequences.append(lines[i])
		if i % 4 == 3:
			labels.append(lines[i])

assert len(sequences) == len(labels)
tot_bind = 0
tot_nonbind = 0

pair_idx = 0
while(pair_idx < len(sequences)):
	i = 0
	for aa in sequences[pair_idx]:
		if labels[pair_idx][i] == "0":
			nonbinding_counts[aa] = nonbinding_counts.get(aa, 0) + 1
			tot_nonbind += 1
		elif labels[pair_idx][i] == "1":
			bind_counts[aa] = bind_counts.get(aa, 0) + 1
			tot_bind += 1
		else:
			assert 0
		i += 1
	pair_idx += 1

bind_counts = dict(sorted(bind_counts.items()))
nonbinding_counts = dict(sorted(nonbinding_counts.items()))

print("Binding sites:")
for aa in bind_counts.keys():
	print("\"" + aa + "\" : " + str(round(bind_counts[aa]/tot_bind, 6)) + ",")
print("----------------------------------\nNon-binding sites:")
for aa in nonbinding_counts.keys():
	print("\"" + aa + "\" : " + str(round(nonbinding_counts[aa]/tot_nonbind, 6)) + ",")
