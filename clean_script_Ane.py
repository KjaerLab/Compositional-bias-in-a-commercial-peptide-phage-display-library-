#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 14:00:23 2022

@author: maltejensen
"""
'''
Load data
'''
import pandas as pd
import numpy as np
import os 
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA, TruncatedSVD
from termcolor import colored
from tqdm import tqdm


    
# file_root = '/home/maltejensen/Documents/CMI/phage_display_data/DKS_004_2102_Competition_assay_t230_R1_AA_Clean.csv'
# file_root = '/home/maltejensen/Documents/CMI/phage_display_data/Ane data/10uL_AA_Clean.csv'
file_root = '/home/maltejensen/Documents/CMI/phage_display_data/Ane data/15_No_X_AA_Clean.csv'

# file_root = '/home/maltejensen/Documents/CMI/phage_display_data/data_2.csv'

df = pd.read_csv(file_root)
# df = pd.read_csv(file_root, delimiter=',', chunksize=10000)
# df = df.get_chunk()


#%%
'''
Onehot encode data before analysis
'''

AA_letters = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

AA_map = {AA_letters[i]: i for i in range(len(AA_letters))}
AA_map_rev = {i: AA_letters[i] for i in range(len(AA_letters))}


# map of properties
prop_map = {
    0: ['G','A','V','L','I','P'],       # Aliphatic side-chains
    1: ['S', 'T'],                      # Polar neutral side-chains
    2: ['N', 'Q'],                      # Amide side-chains
    3: ['C', 'M'],                      # Sulfur-containing side-chains
    4: ['F', 'Y', 'W'],                 # Aromatic side-chains
    5: ['D', 'E'],                      # Anionic side-chains
    6: ['H', 'K', 'R']                  # Cationic side-chains
    }

prop_map_rev = {
    0: 'Aliphatic',
    1: 'Polar',
    2: 'Amide',
    3: 'Sulfur',
    4: 'Aromatic',
    5: 'Anionic',
    6: 'Cationic'
    }

prop_map = {
    0: ['G','A','V','L','I','P'],       # Aliphatic side-chains
    1: ['S', 'T'],                      # Polar neutral side-chains
    2: ['N', 'Q'],                      # Amide side-chains
    3: ['C', 'M'],                      # Sulfur-containing side-chains
    4: ['F', 'Y', 'W'],                 # Aromatic side-chains
    5: ['D', 'E'],                      # Anionic side-chains
    6: ['H', 'K', 'R']                  # Cationic side-chains
    }

prop_map_rev_groups = {
    'Aliphatic': ['G','A','V','L','I','P'],       # Aliphatic side-chains
    'Polar': ['S', 'T'],                      # Polar neutral side-chains
    'Amide': ['N', 'Q'],                      # Amide side-chains
    'Sulfur': ['C', 'M'],                      # Sulfur-containing side-chains
    'Aromatic': ['F', 'Y', 'W'],                 # Aromatic side-chains
    'Anionic': ['D', 'E'],                      # Anionic side-chains
    'Cationic': ['H', 'K', 'R']                  # Cationic side-chains
    }

def OneHotEncodeAA(seq, AA_map, num_seq, type_):
    
    seq = seq.strip('\'')
    crop_seq = seq[:num_seq]
    
    num_col = len(AA_map)
    
    one_hot = np.zeros(num_col*len(crop_seq), dtype=np.int8)
    
    if type_ == 'AA':
        for i in range(len(crop_seq)):
            one_hot[i*num_col + AA_map[crop_seq[i]]] = 1
    elif type_ == 'prop':
        for i in range(len(crop_seq)):
            for key in AA_map.keys():
                if crop_seq[i] in AA_map[key]:
                    one_hot[i*num_col + key] = 1
                    break
    return one_hot


# embedding type
embedding_type = 'AA' # AA or prop

# set number of samples to include and how long the main sequence is
# num_samples = int(1e3)
num_samples = df.shape[0]
num_seq = 12

if embedding_type == 'AA':
    seq_map = AA_map
    seq_map_rev = AA_map_rev
    type_ = 'AA'
elif embedding_type == 'prop':
    seq_map = prop_map
    seq_map_rev = prop_map_rev
    type_ = 'prop'

num_groups = len(seq_map)
one_hot = []

for index, row in tqdm(df.iterrows(), total=len(df)):
    one_hot.append(OneHotEncodeAA(row['AA'], seq_map, num_seq=num_seq, type_=type_))
    
    if index == num_samples-1:
        break

one_hot = np.vstack(one_hot)

y = df['Relative (%)'][:num_samples]

#%%
# import prince

# Perform PCA
pca = PCA(n_components=num_samples if num_samples<one_hot.shape[1] else one_hot.shape[1], svd_solver='full')
# pca = TruncatedSVD(n_components=num_samples if num_samples<one_hot.shape[1] else one_hot.shape[1]-1,
#                     algorithm = 'randomized')

pca.fit(one_hot)    
# pca.fit(one_hot.astype(float)/np.array(list(num_codons.values())*num_seq).reshape(1,-1))
# pca.fit(one_hot.astype(float)/np.sqrt(np.array(list(num_codons.values())*num_seq)).reshape(1,-1))
# pca.fit(one_hot.astype(float)/one_hot.std(axis=0).reshape(1,-1))

pc = pca.transform(one_hot)
# pc = pca.transform(one_hot.astype(float)/np.array(list(num_codons.values())*num_seq).reshape(1,-1) )
# pc = pca.transform(one_hot.astype(float)/np.sqrt(np.array(list(num_codons.values())*num_seq)).reshape(1,-1) )
# pc = pca.transform(one_hot.astype(float)/one_hot.std(axis=0).reshape(1,-1))

#%%
'''
Plot PCA against number of hits (y)
'''
plt.figure()
plt.plot(np.cumsum(pca.explained_variance_ratio_))
plt.grid()
plt.title('explained variance')
plt.show()

# plots that show if component 1 or 2 are correlated in the number of hits (y)
plt.figure()
plt.plot(pc[:,0],y,'.')
plt.grid()
plt.title('PC 1 vs y')
plt.ylim(-0.01, y[1])
plt.xlabel('PC 1')
plt.ylabel('Number of hits')
plt.show()

plt.figure()
plt.plot(pc[:,1],y,'.')
plt.grid()
plt.title('PC 2 vs y')
plt.xlabel('PC 2')
plt.ylabel('Number of hits')
plt.ylim(-0.01, y[1])
plt.show()

#%%
'''
Plots PC 1 vs PC 2 and colors the points after how many hits a particular sequence has
'''

x_pc = 0
y_pc = x_pc+1


ss = int(1e5)

plt.figure(figsize=(15,10))
if pc.shape[0] > ss:
    print('Plotted subsample of data')
    subsample = np.random.choice(pc.shape[0], ss, replace=False)
    subsample = np.sort(subsample)
    plt.scatter(pc[subsample,x_pc][::-1],pc[subsample,y_pc][::-1], s=20, c=np.log(y[subsample][::-1]), cmap='jet')
    # plt.scatter(pc[subsample,0],pc[subsample,1], s=20, edgecolors='black')
else:    
    plt.scatter(pc[:,0][::-1],pc[:,1][::-1], s=20, c=np.log(y[:][::-1]), cmap='jet')


plt.title('PC {} vs PC {}'.format(x_pc+1, y_pc+1))
plt.xlabel('PC ' + str(x_pc+1))
plt.ylabel('PC ' + str(y_pc+1))
plt.grid()
cbar = plt.colorbar(label='Percent')
ticks = cbar.ax.get_yticks()
# points are colored after log(counts), so they need to be converted back
cbar.ax.set_yticklabels(np.round(np.exp(ticks)*100, 3))
plt.show()

#%%
'''
Sum up the occurence of amino acids 
'''

# notice this is for a sequence of 12 (3,4) should be changed for more/less plots
fig, axx = plt.subplots(3,4, figsize=(18,10))

# change this if the y-axis is too small/large
ylim = [0, 0.20]

for position in range(num_seq):
    # axx[position//4,position % 4].bar(list(seq_map_rev.values()), 
    #         one_hot[:,len(seq_map_rev)*position:len(seq_map_rev)*(position+1)].sum(axis=0)/one_hot.shape[0])
    
    axx[position//4,position % 4].bar(list(seq_map_rev.values()), 
            one_hot[:,len(seq_map_rev)*position:len(seq_map_rev)*(position+1)].var(axis=0))
    
    axx[position//4,position % 4].set_title('Position {}'.format(position+1))
    axx[position//4,position % 4].set_ylim(*ylim)
    # baseline
    axx[position//4,position % 4].plot([-1, len(seq_map_rev)+1], [1/len(seq_map_rev), 1/len(seq_map_rev)], '--', color='gray', alpha=0.7)
    axx[position//4,position % 4].set_xlim([-1, len(seq_map_rev)])
    # axx[position//4,position % 4].grid(axis='y')
    
    
    if position % 4 == 0:
        axx[position//4,position % 4].set_ylabel('Fraction')
    
    if embedding_type == 'prop':
        if position<8:
            axx[position//4,position % 4].xaxis.set_ticklabels([])
            axx[position//4,position % 4].yaxis.set_ticklabels([])
            axx[position//4,position % 4].set_xticks([])
            axx[position//4,position % 4].set_yticks([])
        else:
            axx[position//4,position % 4].xaxis.set_tick_params(rotation=-45)
            axx[position//4,position % 4].set_xticklabels(list(seq_map_rev.values()), ha='left', rotation_mode = 'anchor')
    else:
        pass 
plt.show()
#%%
'''
Barplot of global counts of amino acids
'''

AA_count = np.zeros(len(seq_map_rev))

for position in range(num_seq):
    AA_count += one_hot[:,len(seq_map_rev)*position:len(seq_map_rev)*(position+1)].sum(axis=0)
    
AA_count = AA_count/(one_hot.shape[0] * num_seq)

plt.figure()
plt.bar(list(seq_map_rev.values()), AA_count) 
plt.plot([-1, len(seq_map_rev)+1], [1/len(seq_map_rev), 1/len(seq_map_rev)], '--', color='gray', alpha=0.7)
plt.xlim([-1, len(seq_map_rev)])
plt.title('Total amino acid count')
plt.xlabel('Amino acids')
plt.ylabel('Frequency')
plt.show()

#%%
'''
Heatmap of amino acids 
'''
from mpl_toolkits.axes_grid1 import make_axes_locatable

# define baseline
num_codons = {
    'A': 2,
    'R': 3,
    'N': 1,
    'D': 1,
    'C': 1,
    'Q': 2,
    'E': 1,
    'G': 2,
    'H': 1,
    'I': 1,
    'L': 3,
    'K': 1,
    'M': 1,
    'F': 1,
    'P': 2,
    'S': 3,
    'T': 2,
    'W': 1,
    'Y': 1,
    'V': 2
    }

# use these paramters to control if the loading should be in chunks
# and the size of the chunks 
use_chunks = False
chunk_size = 1e5

baseline_freq = np.array([num_codons[aa] for aa in AA_letters]).reshape(-1,1)
baseline_freq = baseline_freq/baseline_freq.sum()

AA_heatmap = []

if use_chunks:
    AA_heatmap = np.zeros((len(seq_map_rev), num_seq))
       
    # Read a chunk
    for i, df_chunk in enumerate(pd.read_csv(file_root, chunksize=chunk_size)):
        if (i + 1) % 5 == 0:
            print('Read {:.1e} rows'.format((i+1)*chunk_size))
    
        one_hot_chunk = []
        # create the one_hot encoding for this chunk
        for row_num, (index, row) in enumerate(df_chunk.iterrows()):
            one_hot_chunk.append(OneHotEncodeAA(row['AA'], seq_map, num_seq=num_seq, type_=type_))   
            
        one_hot_chunk = np.vstack(one_hot_chunk)

        # sum up the heat map
        for position in range(num_seq):
            AA_count = one_hot_chunk[:,len(seq_map_rev)*position:len(seq_map_rev)*(position+1)].sum(axis=0)
            AA_heatmap[:,position] = AA_heatmap[:,position] + AA_count
        
    num_hits = AA_heatmap[:,0].sum()
    AA_heatmap = (AA_heatmap/num_hits) - baseline_freq

else:
    for position in range(num_seq):
        AA_count = one_hot[:,len(seq_map_rev)*position:len(seq_map_rev)*(position+1)].sum(axis=0)/one_hot.shape[0]
        AA_heatmap.append(AA_count.reshape(-1,1))
    
    AA_heatmap = np.hstack(AA_heatmap)-baseline_freq


fig, axx = plt.subplots()
abs_max = np.max(np.abs(AA_heatmap))
img = axx.imshow(AA_heatmap, cmap='bwr', vmin=-abs_max, vmax=abs_max)

axx.set_yticks(np.arange(len(seq_map_rev)))
axx.set_yticklabels(list(seq_map_rev.values()))
axx.set_xticks(np.arange(num_seq))
axx.set_xticklabels(np.arange(1, num_seq+1))

axx.set_xlabel('Position')
axx.set_ylabel('Amino acid')
axx.set_title('Heat map of amino acid position')

divider = make_axes_locatable(axx)
cax = divider.append_axes("right", size="5%", pad=0.1)
cbar = plt.colorbar(img, cax=cax, label='Fraction (difference from baseline)')

plt.show()


#%%
'''
Extract the meaning of the loadings using amino acids
'''
num_loadings = 10
color_cutoff = num_loadings
num_seqs_to_show = 50
pc_to_show = 1 # not zero-indexed
loadings = pca.components_[pc_to_show-1,:]
# loadings = pca.components_[:20,:].mean(axis=0)
load_idx = np.argsort(-loadings)[:num_loadings]

placements = {i: [] for i in range(one_hot.shape[1]//20)}

print('*'*5+' PC ' + str(pc_to_show) +' ' + '*'*5)
print('#'*5 +' Positive influence ' + '#'*5)
print('{:<10}{:<15}{:<10}'.format('Position', 'AA', 'Loading'))

pos_load = {}
for i, l_idx in enumerate(load_idx):
    # find the position
    idx_ = l_idx//num_groups
    AA_ = seq_map_rev[l_idx%num_groups]
    load_ = loadings[l_idx]
    # save the top loadings for coloring output
    if i < color_cutoff:
        pos_load[idx_] = AA_
    
    # placements[idx].append([AA_map_rev[l_idx%20], loadings[l_idx]])
    print('{:<10}{:<15}{:<10.3f}'.format(idx_+1, AA_, load_))

load_idx = np.argsort(loadings)[:num_loadings]

placements = {i: [] for i in range(one_hot.shape[1]//20)}

print()
print('#'*5 +' Negative influence ' + '#'*5)
print('{:<10}{:<15}{:<10}'.format('Position', 'AA', 'Loading'))

neg_load = {}
for i, l_idx in enumerate(load_idx):
    # find the position
    idx_ = l_idx//num_groups
    AA_ = seq_map_rev[l_idx%num_groups]
    load_ = loadings[l_idx]
    # save the top loadings for coloring output
    if i < color_cutoff:
        neg_load[idx_] = AA_
        
    # placements[idx].append([AA_map_rev[l_idx%20], loadings[l_idx]])
    print('{:<10}{:<15}{:<10.3f}'.format(idx_+1, AA_, load_))



print('\n###### Top binders ######')
# print(df.iloc[:20,:])

print('{:<20}{:10}{}'.format(*df.columns))
for i in range(num_seqs_to_show):
    # split AA at every 5 to a space for readability
    # AA_list = [df['AA'][i][j:j+5] for j in np.arange(0,len(df['AA'][i]),5)]
    # print('{:<25}{:<10}{:.4f}'.format(' '.join(AA_list), 
    #                                   df['Number'][i], df['Relative (%)'][i]))
    AA_str = ''
    for idx, AA in enumerate(df['AA'][i][:num_seq]):
        # determine whether any AA sequence matches the top loadings
        if embedding_type == 'AA':
            if idx in pos_load and AA == pos_load[idx]:
                AA_str += colored(AA, color='red', attrs=['bold'])
            elif idx in neg_load and AA == neg_load[idx]:
                AA_str += colored(AA, color='blue', attrs=['bold'])
            else:
                AA_str += AA
        elif embedding_type == 'prop':    
            if idx in pos_load and AA in prop_map_rev_groups[pos_load[idx]]:
                AA_str += colored(AA, color='red', attrs=['bold'])
            elif idx in neg_load and AA in prop_map_rev_groups[neg_load[idx]]:
                AA_str += colored(AA, color='blue', attrs=['bold'])
            else:
                AA_str += AA
    print('{}{}{:<10}{:.4f}'.format(AA_str,' '*(20-num_seq),
                                  df['Number'][i], df['Relative (%)'][i]))
   
#%%
'''
Make fancy histogram with boxes 
'''
from matplotlib.patches import Rectangle
import matplotlib

bins_ = np.append(np.logspace(-5,-2, num=11),1) # Like in the paper
# bins_ = np.append(np.hstack([ [1*10**i,5*10**i] for i in range(3)])*1e-5,
#                   [1e-2, 1])

print(bins_)

binned_data = y.value_counts(bins=bins_)
binned_data_frac = binned_data/len(y)

print(binned_data)
print('')
print(binned_data_frac)


norm = matplotlib.colors.Normalize(vmin=1, 
                                   vmax=len(bins_)-2)

cmap_col = matplotlib.cm.get_cmap('jet')
cumsum_frac = 0
fig, axx = plt.subplots(1)

for i, (count, frac) in enumerate(zip(binned_data, binned_data_frac)):
    if i == 0:
        col = 'gray'
    else:
        col = cmap_col(norm(i))
    
    box = Rectangle((0, cumsum_frac), width=count, height=frac,
                    facecolor = col,
                    edgecolor='black')
    axx.add_patch(box)
    
    cumsum_frac += frac

axx.set_xlim((0, max(binned_data) + max(binned_data)*0.1 ))
axx.set_ylim((0,1))

axx.set_xlabel('Unique seq')
axx.set_ylabel('Fraction of total sequences')

# divider = make_axes_locatable(axx)
# cax = divider.append_axes("right", size="5%", pad=0.1)
# cbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_col), ax=axx,
#               label='Sequence abundance')

# ticks = cbar.ax.get_yticks()


# cbar.ax.set_yticklabels(np.round(np.exp(ticks)*100, 3))


plt.show()
    
#%%
'''
GMM clustering of PCA data
'''
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
import matplotlib as mpl

ss = int(1e5)

x_pc = 0
y_pc = x_pc+1

gmm_components = 6


if pc.shape[0] > ss:
    print('Used subset for estimation')
    subsample = np.random.choice(pc.shape[0], int(ss), replace=False)
    subsample = np.sort(subsample)
    
    gm = GaussianMixture(n_components=gmm_components, random_state=None, n_init=3).fit(pc[subsample,x_pc:y_pc+1])
    gm_means = gm.means_
else:
    gm = GaussianMixture(n_components=gmm_components, random_state=None, n_init=3, verbose=2).fit(pc[:,x_pc:y_pc+1])
    gm_means = gm.means_

pred_ = gm.predict(pc[:,x_pc:y_pc+1])



def make_ellipses(gmm, ax, cmap='Pastel1'):
    cmap_ = mpl.cm.get_cmap(cmap)
    
    for n in range(gmm.n_components):
        if gmm.covariance_type == "full":
            covariances = gmm.covariances_[n][:2, :2]
        elif gmm.covariance_type == "tied":
            covariances = gmm.covariances_[:2, :2]
        elif gmm.covariance_type == "diag":
            covariances = np.diag(gmm.covariances_[n][:2])
        elif gmm.covariance_type == "spherical":
            covariances = np.eye(gmm.means_.shape[1]) * gmm.covariances_[n]
        v, w = np.linalg.eigh(covariances*2)
        u = w[0] / np.linalg.norm(w[0])
        angle = np.arctan2(u[1], u[0])
        angle = 180 * angle / np.pi  # convert to degrees
        v = 2.0 * np.sqrt(2.0) * np.sqrt(v)
        ell = mpl.patches.Ellipse(
            gmm.means_[n, :2], v[0], v[1], 180 + angle, color=cmap_(i)
        )
        ell.set_clip_box(ax.bbox)
        ell.set_alpha(0.5)
        ax.add_artist(ell)
        ax.set_aspect("equal", "datalim")

fig, axx = plt.subplots(1, figsize=(10,7))
axx.grid()


if pc.shape[0] > ss:
    print('Plotted subsample of data')
    subsample = np.random.choice(pc.shape[0], int(ss), replace=False)
    subsample = np.sort(subsample)
    scatter = axx.scatter(pc[subsample,x_pc], pc[subsample,y_pc], c=pred_[subsample], s=20, cmap='tab10')
    axx.legend(handles=scatter.legend_elements()[0], labels=['Cluster {}'.format(i+1) for i in range(gm.n_components)],
               loc='upper right')
    # plt.plot(gm_means[:,0], gm_means[:,1], 'rx', markersize=20, markeredgewidth=4)
   
    for i in range(gm_means.shape[0]):
        axx.text(gm_means[i,0], gm_means[i,1], s=str(i+1))
        
    make_ellipses(gm, axx)
else:    
    scatter = plt.scatter(pc[:,x_pc],pc[:,y_pc], c=pred_, s=20)
    # plt.plot(gm_means[:,0], gm_means[:,1], 'rx', markersize=20, markeredgewidth=4)
    plt.legend(handles=scatter.legend_elements()[0], labels=['Cluster {}'.format(i+1) for i in range(gm.n_components)],
               loc='upper right')



axx.set_title('PC {} vs PC {} clusters'.format(x_pc+1, y_pc+1))
axx.set_xlabel('PC ' + str(x_pc +1))
axx.set_ylabel('PC ' + str(y_pc +1))
plt.show()

#%%
'''
Make heatmaps of the clusters
'''


baseline_freq = np.array([num_codons[aa] for aa in AA_letters]).reshape(-1,1)
baseline_freq = baseline_freq/baseline_freq.sum()

rows = int(gm.n_components**0.5)
cols = int(np.ceil(gm.n_components/int(gm.n_components**0.5)))

fig, axx = plt.subplots(rows, cols,
                        figsize=(10,10))

for c in range(gm.n_components):
    AA_heatmap = []
    
    class_idx = pred_ == c
        
    for position in range(num_seq):
        position_slice = one_hot[class_idx,len(seq_map_rev)*position:len(seq_map_rev)*(position+1)].copy()
        position_slice[position_slice > 0] = 1
        AA_count = position_slice.sum(axis=0)/position_slice.shape[0]
        AA_heatmap.append(AA_count.reshape(-1,1))
    
    AA_heatmap = np.hstack(AA_heatmap)-baseline_freq
    
    if rows == 1:
        axx_idx = (c,)
    else:
        axx_idx = (c//cols, c%cols)
    
    abs_max = np.max(np.abs(AA_heatmap))    
    img = axx[axx_idx].imshow(AA_heatmap, cmap='bwr', vmin=-abs_max, vmax=abs_max)
    
    axx[axx_idx].set_yticks(np.arange(len(seq_map_rev)))
    axx[axx_idx].set_yticklabels(list(seq_map_rev.values()))
    axx[axx_idx].set_xticks(np.arange(num_seq))
    axx[axx_idx].set_xticklabels(np.arange(1, num_seq+1))
    
    axx[axx_idx].set_xlabel('Position')
    axx[axx_idx].set_ylabel('Amino acid')
    axx[axx_idx].set_title('Cluster {} ({:.1e} unique seq, {:.1f}%)'.format(c+1, class_idx.sum(),
                                                                              (class_idx.sum()/one_hot.shape[0])*100 ))
    
    divider = make_axes_locatable(axx[axx_idx])
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = plt.colorbar(img, cax=cax, label='Fraction (difference from baseline)')

plt.tight_layout()    
plt.show()
    
    
#%%n
'''
Simulate data
'''

# define baseline
num_codons = {
    'A': 2,
    'R': 3,
    'N': 1,
    'D': 1,
    'C': 1,
    'Q': 2,
    'E': 1,
    'G': 2,
    'H': 1,
    'I': 1,
    'L': 3,
    'K': 1,
    'M': 1,
    'F': 1,
    'P': 2,
    'S': 3,
    'T': 2,
    'W': 1,
    'Y': 1,
    'V': 2
    }

# use these paramters to control if the loading should be in chunks
# and the size of the chunks 
use_chunks = False
chunk_size = 1e5
num_sim = int(3.7e6)

baseline_freq = np.array([num_codons[aa] for aa in AA_letters])
baseline_freq = baseline_freq/baseline_freq.sum()

one_hot = []

pos_vec = np.array([20*i for i in range(num_seq)])

print('Simulating data..')
for i in tqdm(range(num_sim)): 
    AA_sample = np.random.choice(len(AA_letters), size=num_seq, replace=True, p=baseline_freq)
    oh = np.zeros(len(AA_letters)*num_seq)
    
    oh[AA_sample + pos_vec] = 1
    one_hot.append(oh)
    
one_hot = np.vstack(one_hot)    

y = np.ones(one_hot.shape[0])

#%%
'''
Display confidence intervals for number of read 
'''
import scipy

df['probability'] = df['Relative (%)']/100

num_reads = df['Number'].sum()
# critical value
alpha = 0.05
num_points = 30

critical_z = scipy.stats.norm.ppf(1-(alpha/num_points)/2)
# critical_z = 1.96

ci_reads = critical_z*(((df['probability']*(1-df['probability']))/num_reads)**0.5)

plt.figure()
plt.errorbar(np.arange(1,num_points+1), df['Relative (%)'][:num_points], 
              yerr=ci_reads[:num_points]*100, fmt = '.', barsabove=False, capsize=3, ecolor='gray')
plt.grid()
plt.xticks([1,] + list(np.arange(5,num_points+1,5)) )
plt.title('Prevalence of top binders (95% confidence interval)')
plt.xlabel('# Binder') 
plt.ylabel('% of sample')

for i, (p_, ci) in enumerate(zip(df['Relative (%)'][:num_points], ci_reads[:num_points]), start=1):
    print('{}: {:.1e} +/- {:.1e} %'.format(i, p_, ci*100))
    
plt.show()
#%%

'''
Show bionomial distribution of the subsample form the phage pool
'''hvorforimport scipy

df['probability'] = df['Relative (%)']/100

num_reads = df['Number'].sum()
# critical value
alpha = 0.05
num_points = 30

critical_z = scipy.stats.norm.ppf(1-(alpha/num_points)/2)
# critical_z = 1.96

ci_reads = critical_z*(((df['probability']*(1-df['probability']))/num_reads)**0.5)

plt.figure()
plt.errorbar(np.arange(1,num_points+1), df['Relative (%)'][:num_points], 
              yerr=ci_reads[:num_points]*100, fmt = '.', barsabove=False, capsize=3, ecolor='gray')
plt.grid()
plt.xticks([1,] + list(np.arange(5,num_points+1,5)) )
plt.title('Prevalence of top binders (95% confidence interval)')
plt.xlabel('# Binder') 
plt.ylabel('% of sample')

for i, (p_, ci) in enumerate(zip(df['Relative (%)'][:num_points], ci_reads[:num_points]), start=1):
    print('{}: {:.1e} +/- {:.1e} %'.format(i, p_, ci*100))
    
plt.show()
from scipy.stats import binom

pool = 1e9
pool_rel =  1/pool
seq_num = int(4e6)

use_x_or_more = False


for k in range(5):

    if use_x_or_more:  
        if k == 0:
            print('{} copies: {:.2e}'.format(k, binom.pmf(k, seq_num, pool_rel)))
        else:
            print('{} or more copies: {:.2e}'.format(k, 1-binom.cdf(k-1, seq_num, pool_rel)))
    else:
        print('{} copies: {:.5e}'.format(k, binom.pmf(k, seq_num, pool_rel)))


#%%

import umap
import umap.plot
from sklearn.manifold import TSNE
import time
import numba
import numpy as np

# chebyshev is interesting
# mahalanobis
# wminkowski (very slow)

@numba.njit()
def myDist(a, b):
    dot = np.dot(a, b)
    if dot < 2.0:
        dot = 0.0
    
    return 1-dot/12.0

t0 = time.time()
# reducer = umap.UMAP(n_components=2, n_neighbors=30, min_dist=0, metric='cosine', verbose=True)
reducer = umap.UMAP(n_components=2, n_neighbors=10, min_dist=0, metric=myDist, verbose=True)

mapper = reducer.fit(one_hot[:int(1e6), :])
umap_fit = mapper.transform(one_hot[:int(1e6), :])


# mapper = reducer.fit(one_hot.astype(float)/np.sqrt(np.array(list(num_codons.values())*num_seq)).reshape(1,-1) )
# umap_fit = mapper.transform(one_hot.astype(float)/np.sqrt(np.array(list(num_codons.values())*num_seq)).reshape(1,-1) )
# umap_fit = TSNE(n_components=2, learning_rate='auto',
#                                   init='random').fit_transform(one_hot)

print('Time: {:.1f} min.'.format((time.time()-t0)/60 ))

#%%
ss = int(1e5)

plt.figure(figsize=(15,10))
if umap_fit.shape[0] > ss:
    print('Plotted subsample of data')
    subsample = np.random.choice(umap_fit.shape[0], ss, replace=False)
    subsample = np.sort(subsample)
    plt.scatter(umap_fit[subsample,0][::-1],umap_fit[subsample,1][::-1], s=2, c=np.log(y[subsample][::-1]), cmap='jet', alpha=0.7)
    # plt.scatter(pc[subsample,0],pc[subsample,1], s=20, edgecolors='black')
else:    
    plt.scatter(umap_fit[:,0][::-1],umap_fit[:,1][::-1], s=20, c=np.log(y[:][::-1]), cmap='jet')


plt.title('UMAP')
plt.xlabel('axis 1')
plt.ylabel('axis 2')
plt.grid()
cbar = plt.colorbar(label='Percent')
ticks = cbar.ax.get_yticks()
# points are colored after log(counts), so they need to be converted back
cbar.ax.set_yticklabels(np.round(np.exp(ticks)*100, 3))
plt.show()

#%%


# p = umap.plot.interactive(mapper, values=np.log(y), hover_data=df.iloc[:one_hot.shape[0], :], 
#                           point_size=4, cmap='jet', height=1000, width=1400)
# umap.plot.show(p)

#%%
# umap.plot.points(mapper, values=np.log(y), cmap='jet')
umap.plot.points(mapper)
plt.show()

#%%

umap.plot.connectivity(fit_, edge_bundling='hammer')
# umap.plot.connectivity(fit_)
#%%

'''
GMM clustering of UMAP data
'''
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
import matplotlib as mpl

# ss = int(1e5)
ss = int(1e6)

gmm_components = 8 

if umap_fit.shape[0] > ss:
    print('Used subset for estimation')
    subsample = np.random.choice(umap_fit.shape[0], int(ss), replace=False)
    subsample = np.sort(subsample)
    
    gm = GaussianMixture(n_components=gmm_components, random_state=None, n_init=3).fit(umap_fit[subsample,:])
    gm_means = gm.means_
else:
    gm = GaussianMixture(n_components=gmm_components, random_state=None, n_init=1, verbose=2).fit(umap_fit)
    gm_means = gm.means_

pred_ = gm.predict(umap_fit)



def make_ellipses(gmm, ax, cmap='Pastel1'):
    cmap_ = mpl.cm.get_cmap(cmap)
    
    for n in range(gmm.n_components):
        if gmm.covariance_type == "full":
            covariances = gmm.covariances_[n][:2, :2]
        elif gmm.covariance_type == "tied":
            covariances = gmm.covariances_[:2, :2]
        elif gmm.covariance_type == "diag": 
            covariances = np.diag(gmm.covariances_[n][:2])
        elif gmm.covariance_type == "spherical":
            covariances = np.eye(gmm.means_.shape[1]) * gmm.covariances_[n]
        v, w = np.linalg.eigh(covariances*2)
        u = w[0] / np.linalg.norm(w[0])
        angle = np.arctan2(u[1], u[0])
        angle = 180 * angle / np.pi  # convert to degrees
        v = 2.0 * np.sqrt(2.0) * np.sqrt(v)
        ell = mpl.patches.Ellipse(
            gmm.means_[n, :2], v[0], v[1], 180 + angle, color=cmap_(i)
        )
        ell.set_clip_box(ax.bbox)
        ell.set_alpha(0.5)
        ax.add_artist(ell)
        ax.set_aspect("equal", "datalim")



fig, axx = plt.subplots(1, figsize=(10,7))
axx.grid()


if umap_fit.shape[0] > 1e5:
    print('Plotted subsample of data')
    subsample = np.random.choice(umap_fit.shape[0], int(1e5), replace=False)
    subsample = np.sort(subsample)
    scatter = axx.scatter(umap_fit[subsample,0], umap_fit[subsample,1], c=pred_[subsample], s=0.1, cmap='tab10', alpha=1)
    axx.legend(handles=scatter.legend_elements()[0], 
               labels=['Cluster {}'.format(i+1) for i in range(gm.n_components)], loc='upper right')
    # plt.plot(gm_means[:,0], gm_means[:,1], 'rx', markersize=20, markeredgewidth=4)
   
    for i in range(gm_means.shape[0]):
        axx.text(gm_means[i,0], gm_means[i,1], s=str(i+1))
        
    make_ellipses(gm, axx)
else:    
    scatter = plt.scatter(umap_fit[:,0],umap_fit[:,1], c=pred_, s=2)
    # plt.plot(gm_means[:,0], gm_means[:,1], 'rx', markersize=20, markeredgewidth=4)
    plt.legend(handles=scatter.legend_elements()[0], labels=['Cluster {}'.format(i+1) for i in range(gm.n_components)],
               loc='upper right')



# axx.set_title('PC {} vs PC {} clusters'.format(x_pc+1, y_pc+1))
# axx.set_xlabel('PC ' + str(x_pc +1))
# axx.set_ylabel('PC ' + str(y_pc +1))
plt.show()

#%%
'''
Make heatmaps of the clusters
'''
from mpl_toolkits.axes_grid1 import make_axes_locatable

baseline_freq = np.array([num_codons[aa] for aa in AA_letters]).reshape(-1,1)
baseline_freq = baseline_freq/baseline_freq.sum()

rows = int(gm.n_components**0.5)
cols = int(np.ceil(gm.n_components/int(gm.n_components**0.5)))

fig, axx = plt.subplots(rows, cols,
                        figsize=(10,10))
one_hot_temp = one_hot[:int(1e6), :]
for c in range(gm.n_components):
    AA_heatmap = []
    
    class_idx = pred_ == c
        
    for position in range(num_seq):
        position_slice = one_hot_temp[class_idx,len(seq_map_rev)*position:len(seq_map_rev)*(position+1)].copy()
        position_slice[position_slice > 0] = 1
        AA_count = position_slice.sum(axis=0)/position_slice.shape[0]
        AA_heatmap.append(AA_count.reshape(-1,1))
    
    AA_heatmap = np.hstack(AA_heatmap)-baseline_freq
    
    if rows == 1:
        axx_idx = (c,)
    else:
        axx_idx = (c//cols, c%cols)
    
    abs_max = np.max(np.abs(AA_heatmap))    
    img = axx[axx_idx].imshow(AA_heatmap, cmap='bwr', vmin=-abs_max, vmax=abs_max)
    
    axx[axx_idx].set_yticks(np.arange(len(seq_map_rev)))
    axx[axx_idx].set_yticklabels(list(seq_map_rev.values()), fontsize=8)
    axx[axx_idx].set_xticks(np.arange(num_seq))
    axx[axx_idx].set_xticklabels(np.arange(1, num_seq+1),fontsize=8)
    
    axx[axx_idx].set_xlabel('Position', fontsize=8)
    axx[axx_idx].set_ylabel('Amino acid', fontsize=8)
    axx[axx_idx].set_title('Cluster {} ({:.1e}, {:.1f}%)'.format(c+1, class_idx.sum(),
                                                                              (class_idx.sum()/one_hot.shape[0])*100 ),
                           fontsize=8)
    
    divider = make_axes_locatable(axx[axx_idx])
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = plt.colorbar(img, cax=cax, label='Fraction (diff baseline)')

plt.tight_layout()    
plt.show()

#%%
'''
Make PCA on subsets
'''
cluster_subset = one_hot[pred_ == 10-1, :]


#%%
pca_cluster = PCA(n_components=num_samples if num_samples<one_hot.shape[1] else one_hot.shape[1], svd_solver='full')

pca_cluster.fit(cluster_subset)    

pc_subset = pca_cluster.transform(cluster_subset)
#%%

plt.figure()
plt.scatter(pc_subset[:,0], pc_subset[:,1], s=0.1)
plt.show()
#%%
'''
Print PCA for subset
'''
num_loadings = 10
color_cutoff = num_loadings
num_seqs_to_show = 10
pc_to_show = 1 # not zero-indexed
loadings = pca_cluster.components_[pc_to_show-1,:]
# loadings = pca.components_[:20,:].mean(axis=0)
load_idx = np.argsort(-loadings)[:num_loadings]

placements = {i: [] for i in range(one_hot.shape[1]//20)}

print('*'*5+' PC ' + str(pc_to_show) +' ' + '*'*5)
print('#'*5 +' Positive influence ' + '#'*5)
print('{:<10}{:<15}{:<10}'.format('Position', 'AA', 'Loading'))

pos_load = {}
for i, l_idx in enumerate(load_idx):
    # find the position
    idx_ = l_idx//num_groups
    AA_ = seq_map_rev[l_idx%num_groups]
    load_ = loadings[l_idx]
    # save the top loadings for coloring output
    if i < color_cutoff:
        pos_load[idx_] = AA_
    
    # placements[idx].append([AA_map_rev[l_idx%20], loadings[l_idx]])
    print('{:<10}{:<15}{:<10.3f}'.format(idx_+1, AA_, load_))

load_idx = np.argsort(loadings)[:num_loadings]

placements = {i: [] for i in range(one_hot.shape[1]//20)}

print()
print('#'*5 +' Negative influence ' + '#'*5)
print('{:<10}{:<15}{:<10}'.format('Position', 'AA', 'Loading'))

neg_load = {}
for i, l_idx in enumerate(load_idx):
    # find the position
    idx_ = l_idx//num_groups
    AA_ = seq_map_rev[l_idx%num_groups]
    load_ = loadings[l_idx]
    # save the top loadings for coloring output
    if i < color_cutoff:
        neg_load[idx_] = AA_
        
    # placements[idx].append([AA_map_rev[l_idx%20], loadings[l_idx]])
    print('{:<10}{:<15}{:<10.3f}'.format(idx_+1, AA_, load_))



#%%
recon = mapper.inverse_transform(np.array([4.3, 8.7]).reshape(1,-1))
# recon = fit_.inverse_transform(np.array([12.7,]).reshape(1,-1))


AA_heatmap = []


for position in range(num_seq):
    AA_count = recon[0,len(seq_map_rev)*position:len(seq_map_rev)*(position+1)]
    AA_heatmap.append(AA_count.reshape(-1,1))
    # AA_count = one_hot[:,len(seq_map_rev)*position:len(seq_map_rev)*(position+1)].sum(axis=0)/one_hot.shape[0]
    # AA_heatmap.append(AA_count.reshape(-1,1))

AA_heatmap = np.hstack(AA_heatmap)

fig, axx = plt.subplots(figsize=(5,8))
# abs_max = np.max(np.abs(AA_heatmap))
img = axx.imshow(AA_heatmap, cmap='viridis')
# img = axx.imshow(AA_heatmap, cmap='bwr', vmin=-abs_max, vmax=abs_max)

axx.set_yticks(np.arange(len(seq_map_rev)))
axx.set_yticklabels(list(seq_map_rev.values()))
axx.set_xticks(np.arange(num_seq))
axx.set_xticklabels(np.arange(1, num_seq+1))

axx.set_xlabel('Position')
axx.set_ylabel('Amino acid')
axx.set_title('Heat map of amino acid position')

divider = make_axes_locatable(axx)
cax = divider.append_axes("right", size="5%", pad=0.1)
cbar = plt.colorbar(img, cax=cax, label='Fraction (difference from baseline)')

plt.show()



#%%
'''
Plots PC 1 vs PC 2 and colors the points after how many hits a particular sequence has
'''

x_pc = 0
y_pc = 1


ss = int(1e5)

plt.figure(figsize=(15,10))
if umap_fit.shape[0] > ss:
    print('Plotted subsample of data')
    subsample = np.random.choice(pc.shape[0], ss, replace=False)
    subsample = np.sort(subsample)
    # plt.scatter(pc[subsample,x_pc][::-1],pc[subsample,y_pc][::-1], s=20, c=np.log(y[subsample][::-1]), cmap='jet')
    plt.scatter(umap_fit[subsample,0],umap_fit[subsample,1], s=20, c=np.log(y[subsample][::-1]), cmap='jet')
    
    # plt.scatter(pc[subsample,0],pc[subsample,1], s=20, edgecolors='black')
else:    
    plt.scatter(umap_fit[subsample,0],umap_fit[subsample,1], s=20, c=np.log(y[:][::-1]), cmap='jet')


plt.title('PC {} vs PC {}'.format(x_pc+1, y_pc+1))
plt.xlabel('PC ' + str(x_pc+1))
plt.ylabel('PC ' + str(y_pc+1))
plt.grid()
cbar = plt.colorbar(label='Percent')2
ticks = cbar.ax.get_yticks()
# points are colored after log(counts), so they need to be converted back
cbar.ax.set_yticklabels(np.round(np.exp(ticks)*100, 3))
plt.show()


#%%
'''
Load several rounds of biopanning and compare supernatant with the eluant
'''
import os
import pandas as pd
import numpy as np
from tqdm import tqdm

data_root = '/mnt/SUND/Study database/Anders Wilgaard Sinkjær (AWS)/AWS_005_2206_Ph.D.12 hCD4/AWS_005_2206_Ph.D.12 hCD4 Round 1/Analyzed data'

sample_1 = '/mnt/SUND/Study database/Ane Beth Sloth (ABS)/ABS_011_2206_SAPhD12.v2/Analyzed data/SA/8_SA_AA_Clean.csv'
sample_2 = '/mnt/SUND/Study database/Ane Beth Sloth (ABS)/ABS_011_2206_SAPhD12.v2/Analyzed data/SA/1_SA_AA_Clean.csv'


df_1 = pd.read_csv(sample_1)
df_2 = pd.read_csv(sample_2)


''' Load data '''

elu_num = 50000
sup_num = 100

elu_round_1_df = df_1.iloc[:elu_num,:]
elu_round_2_df = df_2.iloc[:elu_num,:]


sup_round_1_df = pd.read_csv(sup_rep_1_round_1_root, nrows=sup_num)
sup_round_2_df = pd.read_csv(sup_rep_1_round_2_root, nrows=sup_num)



''' sort data by rules '''
# necessary increase
elu_incr = 1.0
use_to_100_sup = False

results = {
        'AA': [],
        'round 1 number': [],
        'round 2 number': [],
        'round 1 %': [],
        'round 2 %': [],
        'round 1 index': [],
        'round 2 index': [],
    }

not_incr_count = 0
not_in_top = 0
in_sup_count = 0

for idx, row in tqdm(elu_round_1_df.iterrows()):
    # check if the round 1 target is in the round 2
    tmp_idx_elu = np.where(elu_round_2_df['AA'] == row['AA'])[0]
    
    if len(tmp_idx_elu) > 1:
        raise ValueError('Found more than one copy in round 2\nat index: ' + str(tmp_idx_elu))
    elif len(tmp_idx_elu) == 0:
        not_in_top += 1
    else:
        # row round 2
        row_2 = elu_round_2_df.loc[tmp_idx_elu[0]]
        # check if the found target in round 2 is also in the sup round 2
        sup_idx_2 = np.where(sup_round_2_df['AA'] == row_2['AA'])[0]
        
        
        # rel_inc_bool = row_2['Relative (%)']/row['Relative (%)'] >= elu_incr
        rel_inc_bool = True
        
        if use_to_100_sup:
            in_sup_2_bool = len(sup_idx_2) == 0
        else:
            in_sup_2_bool = True
        
        if  rel_inc_bool and in_sup_2_bool:
            results['AA'].append(row['AA'])
            
            results['round 1 number'].append(row['Number'])
            results['round 2 number'].append(row_2['Number'])
            
            results['round 1 %'].append(row['Relative (%)'])
            results['round 2 %'].append(row_2['Relative (%)'])
            
            results['round 1 index'].append(idx+1)
            results['round 2 index'].append(tmp_idx_elu[0]+1)
        else:
            not_incr_count += 0 if rel_inc_bool else 1
            in_sup_count += 0 if in_sup_2_bool else 1
            
        
results_df = pd.DataFrame(results)
results_df = results_df.sort_values(by='round 2 %', ascending=False).reset_index(drop=True)
results_df['rel diff'] = results_df['round 1 %']/results_df['round 2 %']

print('')
print('Number of targets kept: {}'.format(results_df.shape[0]))
print('')
print('Number not in top {}: {}'.format(elu_num, not_incr_count))
print('Number below threshold: {}'.format(not_incr_count))
print('Number being in the supernatant: {}'.format(in_sup_count))

#%%
plt.plot(1/results_df['rel diff'],'.')
N = 100
mov_avg = np.convolve(1/results_df['rel diff'], np.ones(N)/N, mode='valid')
plt.plot(mov_avg)
plt.show()

#%%
plt.figure()
plt.plot(results_df['round 1 %'], results_df['round 2 %'],'.')
plt.plot([0,2],[0,2])
plt.show()

#%%
high_ampl_2 = results_df[results_df['rel diff'] < 0.5]
high_ampl_1 = result_2[result_2['rel diff'] < 0.5]

#%%
import pandas as pd
# sample_1 = '/mnt/SUND/Study database/Ane Beth Sloth (ABS)/ABS_011_2206_SAPhD12.v2/Analyzed data/SA/1_SA_AA_Clean.csv'
# sample_2 = '/mnt/SUND/Study database/Ane Beth Sloth (ABS)/ABS_011_2206_SAPhD12.v2/Analyzed data/SA/2_SA_AA_Clean.csv'
# sample_1 = '/mnt/SUND/Study database/Ane Beth Sloth (ABS)/ABS_011_2206_SAPhD12.v2/Analyzed data/SA/5_SA_AA_Clean.csv'
# sample_2 = '/mnt/SUND/Study database/Ane Beth Sloth (ABS)/ABS_011_2206_SAPhD12.v2/Analyzed data/SA/6_SA_AA_Clean.csv'


#%%
import scipy.stats
import matplotlib.pyplot as plt

df= df_1

df['probability'] = df['Relative (%)']/100

num_reads = df['Number'].sum()
# critical value
alpha = 0.05
num_points = 30

critical_z = scipy.stats.norm.ppf(1-(alpha/num_points)/2)
# critical_z = 1.96

ci_reads = critical_z*(((df['probability']*(1-df['probability']))/num_reads)**0.5)

plt.figure()
plt.errorbar(np.arange(1,num_points+1), df['Relative (%)'][:num_points], 
              yerr=ci_reads[:num_points]*100, fmt = '.', barsabove=False, capsize=3, ecolor='gray')
plt.grid(color='gray', axis='y')
plt.xticks([1,] + list(np.arange(5,num_points+1,5)) )
plt.title('Prevalence of top binders (95% confidence interval)')
plt.xlabel('# Binder') 
plt.ylabel('% of sample')

for i, (p_, ci) in enumerate(zip(df['Relative (%)'][:num_points], ci_reads[:num_points]), start=1):
    print('{}: {:.1e} +/- {:.1e} %'.format(i, p_, ci*100))
    
plt.show()

plt.figure()
plt.plot( (ci_reads[:num_points]*100)/df['Relative (%)'][:num_points])
#%%

out = '/home/maltejensen/Documents/CMI/phage_display_data/Ane data/SAPphD12v2'

file_path = '/home/maltejensen/Documents/CMI/phage_display_data/Ane data/SAPphD12v2/4_SA_AA_Clean.csv'

new_file = []

with open(file_path, 'r') as f:
    for line in f:
        new_file.append(line.split(';')[0])

with open(out + '/4_SA_AA_Clean_new.csv', 'w') as f:
    for line in new_file:
        f.write(line+'\n')

#%%
for idx, row in tqdm(df_1.iterrows()):
    if len(np.where(row['AA'] == df_2['AA'])[0]) != 0:
        print(idx, row)




