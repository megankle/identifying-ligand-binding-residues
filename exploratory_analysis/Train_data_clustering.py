#!/usr/bin/env python
# coding: utf-8

# In[1]:


sequences = []
labels = []


with open('/Users/dharit/Downloads/Bioinfo_clustering_project.txt') as f:
   lines = [line.rstrip() for line in f]
   for i in range(len(lines)):
     if i % 4 == 2:
       sequences.append(lines[i])
     if i % 4 == 3:
       labels.append(lines[i])


# In[2]:


# Finds out the length of each protein in the sequences[]:
sequences = []
labels = []

with open('/Users/dharit/Downloads/Bioinfo_clustering_project.txt') as f:
   lines = [line.rstrip() for line in f]
   for i in range(len(lines)):
     if i % 4 == 2:
       sequences.append(lines[i])
     if i % 4 == 3:
       labels.append(lines[i])

# Get length of each protein sequence
for i in range(3): # for the first 3 sequences
    seq_len = len(sequences[i])
    print(f"The length of sequence {i+1} is {seq_len}")


# In[3]:


# Calculate the digram frequencies for each protein sequence and stores as the list.
sequences = []
labels = []


with open('/Users/dharit/Downloads/Bioinfo_clustering_project.txt') as f:
   lines = [line.rstrip() for line in f]
   for i in range(len(lines)):
     if i % 4 == 2:
       sequences.append(lines[i])
     if i % 4 == 3:
       labels.append(lines[i])
def digram_freq(seq):
    freq = {}
    for i in range(len(seq)-1):
        pair = seq[i:i+2]
        if pair in freq:
            freq[pair] += 1
        else:
            freq[pair] = 1
    # Calculate frequency as percentage
    total_pairs = len(seq) - 1
    for pair in freq:
        freq[pair] = round((freq[pair] / total_pairs) * 100, 3)
    return freq

# Calculate digram frequencies for each sequence and store in a list
freq_list = []
for seq in sequences:
    freq = digram_freq(seq)
    freq_list.append(freq)


# In[4]:


# Calculates the ratio of 1s to total count:
sequences = []
labels = []


with open('/Users/dharit/Downloads/Bioinfo_clustering_project.txt') as f:
   lines = [line.rstrip() for line in f]
   for i in range(len(lines)):
     if i % 4 == 2:
       sequences.append(lines[i])
     if i % 4 == 3:
       labels.append(lines[i])
    
label_freq = {}
for label in labels:
    label_freq[label] = {'0': 0, '1': 0}
    for l in label:
        if l == '0':
            label_freq[label]['0'] += 1
        elif l == '1':
            label_freq[label]['1'] += 1
ratios = []

for label in labels:
    count_0 = label.count('0')
    count_1 = label.count('1')
    total_count = count_0 + count_1
    if total_count > 0:
        ratio = count_1 / total_count
        ratios.append(round(ratio, 3))
    else:
        ratios.append(None)


# In[5]:


import pandas as pd

# Load sequences and labels
sequences = []
labels = []
with open('/Users/dharit/Downloads/Bioinfo_clustering_project.txt') as f:
    lines = [line.rstrip() for line in f]
    for i in range(len(lines)):
        if i % 4 == 2:
            sequences.append(lines[i])
        if i % 4 == 3:
            labels.append(lines[i])

# Get length of each protein sequence
seq_lengths = [len(seq) for seq in sequences]

# Calculate digram frequencies for each sequence and store in a list of dicts
freq_list = []
for seq in sequences:
    freq = digram_freq(seq)
    freq_list.append(freq)

# Create a list of dicts to store the data for each sequence
data = []
for i, seq in enumerate(sequences):
    row = {'Sequence': f'Sequence-{i+1}', 'Length': seq_lengths[i]}
    row.update(freq_list[i])
    label_counts = label_freq[labels[i]]
    row['Proportion of binding sites'] = ratios[i]
    data.append(row)

# Create the data frame
df = pd.DataFrame(data)

# Reorder the columns
cols = ['Sequence', 'Length'] + sorted(list(df.columns[2:-3]))
cols += list(df.columns[-3:])
cols += list(df.columns[-2:-1])
cols += list(df.columns[-1:])
df = df[cols]

# Print the data frame
df.fillna(0, inplace=True)
df.head()


# In[6]:


df_new=df.iloc[:,1:]
df_new.head()
df_new_arr=df_new.values


# In[7]:


df_new_arr


# In[8]:


from sklearn.cluster import KMeans

# Assuming your dataframe is named "df"
# Create an instance of the KMeans class with the desired number of clusters
kmeans = KMeans(n_clusters=14)

# Fit the model to your data
kmeans.fit(df_new_arr)

# Get the predicted labels for your data
labels = kmeans.predict(df_new_arr)

# Get the cluster centers
centers = kmeans.cluster_centers_

print(centers)


# In[18]:


from matplotlib import pyplot as plt

df_new['Clusters'] = kmeans.labels_

plt.scatter(df_new["Length"], df_new["Proportion of binding sites"], c=df_new["Clusters"])
plt.xlabel("Length of protein sequences")
plt.ylabel("Proportion of binding sites")
plt.title("K-means Clustering Results on Train data")
plt.savefig("cluster_plot.png")


# In[10]:


plt.savefig("cluster_plot.png")


# In[11]:


df.fillna(0, inplace=True)
df.head()


# In[12]:


import pandas as pd

# Calculate the digram frequencies for each protein sequence and store as a list
sequences = []
labels = []

with open('/Users/dharit/Downloads/Bioinfo_clustering_project.txt') as f:
   lines = [line.rstrip() for line in f]
   for i in range(len(lines)):
     if i % 4 == 2:
       sequences.append(lines[i])
     if i % 4 == 3:
       labels.append(lines[i])

def digram_freq(seq):
    freq = {}
    for i in range(len(seq)-1):
        pair = seq[i:i+2]
        if pair in freq:
            freq[pair] += 1
        else:
            freq[pair] = 1
    # Calculate frequency as percentage
    total_pairs = len(seq) - 1
    for pair in freq:
        freq[pair] = round((freq[pair] / total_pairs) * 100, 3)
    return freq

# Calculate digram frequencies for each sequence and store in a list
freq_list = []
for seq in sequences:
    freq = digram_freq(seq)
    freq_list.append(freq)

# Sort the digram frequencies in alphabetical order
digrams = sorted(list(set([pair for freq in freq_list for pair in freq])))

# Create a dictionary of the data
data = {
    'Sequence ID': [f'Sequence-{i+1}' for i in range(len(sequences))],
}

for digram in digrams:
    data[digram] = []
    for freq in freq_list:
        if digram in freq:
            data[digram].append(freq[digram])
        else:
            data[digram].append(0)

# Create the dataframe
df_2 = pd.DataFrame(data)

# Print the dataframe
df_2.head()


# In[13]:


df_2.to_csv('df2_new.txt', sep='\t', index=False)


# In[ ]:




