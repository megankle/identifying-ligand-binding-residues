#!/usr/bin/env python3
import numpy as np
import sys
from sklearn import metrics
import matplotlib.pyplot as plt

train_type = sys.argv[1]

sequences = []
labels = []
confusion_matrix = np.zeros((2, 2))

with open('Val.dat') as f:
	lines = [line.rstrip() for line in f]
	for i in range(len(lines)):
		if i % 4 == 2:
			sequences.append(lines[i])
		if i % 4 == 3:
			labels.append(lines[i])

if train_type == 'A':
	# AUGMENTED DATA
	b_probs = {"A" : 0.039503, "C" : 0.091378, "D" : 0.108318, "E" : 0.068246, "F" : 0.047803, "G" : 0.069111, "H" : 0.102306, "I" : 0.030924, "K" : 0.040838, "L" : 0.062511, "M" : 0.020998, "N" : 0.044262, "P" : 0.021394, "Q" : 0.019814, "R" : 0.046551, "S" : 0.045606, "T" : 0.050425, "V" : 0.040416, "W" : 0.016672, "Y" : 0.032925}

	n_probs = {"A" : 0.08322, "C" : 0.011925, "D" : 0.057587, "E" : 0.066587, "F" : 0.040925, "G" : 0.076105, "H" : 0.023373, "I" : 0.057382, "K" : 0.05654, "L" : 0.092476, "M" : 0.022539, "N" : 0.041817, "P" : 0.048513, "Q" : 0.036405, "R" : 0.053327, "S" : 0.05714, "T" : 0.054073, "V" : 0.071185, "W" : 0.014536, "Y" : 0.034346}

	transition_probs = {'00': 0.977962, '01': 0.022038, '11': 0.285805, '10': 0.714195}
	start_probs = {'1': 0.009095, '0': 0.990905}
elif train_type == 'O':
	# ORIGINAL DATA
	b_probs = {"A" : 0.04062, "C" : 0.100625, "D" : 0.127618, "E" : 0.076496, "F" : 0.029424, "G" : 0.06961, "H" : 0.107453, "I" : 0.032172, "K" : 0.041951, "L" : 0.040244, "M" : 0.018487, "N" : 0.045944, "P" : 0.022654, "Q" : 0.021757, "R" : 0.047998, "S" : 0.050862, "T" : 0.051123, "V" : 0.034892, "W" : 0.009866, "Y" : 0.030205}

	n_probs = {"A" : 0.082191, "C" : 0.011537, "D" : 0.05794, "E" : 0.067916, "F" : 0.040825, "G" : 0.074159, "H" : 0.022315, "I" : 0.059186, "K" : 0.057913, "L" : 0.094056, "M" : 0.022264, "N" : 0.042387, "P" : 0.046683, "Q" : 0.03601, "R" : 0.051706, "S" : 0.058774, "T" : 0.053668, "V" : 0.071317, "W" : 0.01382, "Y" : 0.035334}
	transition_probs = {'00': 0.980037, '01': 0.019963, '10': 0.717781, '11': 0.282219}
	start_probs = {'0': 0.988435, '1': 0.011565}
else:
	assert 0

def viterbi(sequence):
	# algorithm explanation: https://www.cis.upenn.edu/~cis2620/notes/Example-Viterbi-DNA.pdf
	vit_matrix = np.zeros((2, len(sequence)))

	vit_matrix[0][0] = np.log2(start_probs['0']) + np.log2(n_probs[sequence[0]])
	vit_matrix[1][0] = np.log2(start_probs['1']) + np.log2(b_probs[sequence[0]])

	for i in range(1, len(sequence)):
		aa = sequence[i]
		vit_matrix[0][i] = np.log2(n_probs[aa]) + max(vit_matrix[0][i-1] + np.log2(transition_probs['00']), vit_matrix[1][i-1] + np.log2(transition_probs['10']))
		vit_matrix[1][i] = np.log2(b_probs[aa]) + max(vit_matrix[0][i-1] + np.log2(transition_probs['01']), vit_matrix[1][i-1] + np.log2(transition_probs['11']))

	# backtrack to recreate answer

	i = len(sequence) - 1

	if(vit_matrix[0][i] > vit_matrix[1][i]):
	   answer = '0'
	else:
	   answer = '1'

	while i > 0:
		# get next state you transitioned to
		temp = '0' if (vit_matrix[0][i-1] + np.log2(transition_probs['0' + answer[-1]])) > (vit_matrix[1][i-1] + np.log2(transition_probs['1' + answer[-1]])) else '1'
		answer += temp
		i -= 1 

	return answer[::-1]

actual = []
predicted = []
i = 0
for seq in sequences:
	if i % 5000 == 0:
		print(i)
	ans = viterbi(seq)
	j = 0
	while(j < len(ans)):
		confusion_matrix[int(labels[i][j])][int(ans[j])] += 1
		actual.append(labels[i][j])
		predicted.append(ans[j])
		j += 1
	i += 1

print(confusion_matrix)

confusion_matrix = metrics.confusion_matrix(actual, predicted)

cm_display = metrics.ConfusionMatrixDisplay(confusion_matrix = confusion_matrix, display_labels = ["Non-binding", "Binding"])

cm_display.plot()
plt.savefig("confusion_matrix.png", bbox_inches="tight")

'''
temp = ""
for let in sequence:
	temp += let + " & "
print(temp)
for i in range(2):
	to_print = "S" if i == 0 else "T"
	for j in range(len(sequence)):
		to_print += " & " + str(round(vit_matrix[i][j], 2))
'''
