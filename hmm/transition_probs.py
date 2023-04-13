#!/usr/bin/env python3

labels = []

counts = dict()
start_counts = dict()

with open('Aug-Train.dat') as f:
	lines = [line.rstrip() for line in f]
	for i in range(len(lines)):
		if i % 4 == 3:
			labels.append(lines[i])

for label in labels:
	i = 0
	while(i < len(label) - 1):
		di = label[i:i+2]
		counts[di] = counts.get(di, 0) + 1
		i += 1
	start_counts[label[0]] = start_counts.get(label[0], 0) + 1

counts = dict(sorted(counts.items()))
print(counts)
print("-------------------------")
print(start_counts)
