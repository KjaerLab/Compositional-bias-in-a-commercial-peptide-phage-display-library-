#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 11:50:43 2022

@author: maltejensen
"""

'''
Parse Fastq data to nucleotides
'''
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import re 
from termcolor import colored
import time

G_codons = ['GGG', 'GGT', 'GGA', 'GGC']
S_codons = ['TCA', 'TCC', 'TCG', 'TCT', 'ACG', 'AGT']

GGGS_seqs = []
for c1 in range(len(G_codons)):
    for c2 in range(len(G_codons)):
        for c3 in range(len(G_codons)):
            for c4 in range(len(S_codons)):
                # print(c1,c2,c3,c4)
                GGGS_seqs.append(G_codons[c1]+G_codons[c2]+G_codons[c3]+S_codons[c4])

                

def hasGGGSEnd(seq, GGGS_seqs):
    return not seq_tmp[-12:] in GGGS_seqs
       
def hasCodon(seq, target_codon):
    codons = [seq[i:i+3] == target_codon for i in range(0,36,3)]
    return any(codons), np.where(codons)[0] + 1 

# def ascii2phred(seq):
#     return np.array([ord(a)-33 for a in seq])

def ascii2phred(seq):
    return np.array(list(map(ord, seq)))-33

def findAll(substr, string):
    return [m.start() for m in re.finditer(substr, string)]

def sortDictByKey(dict_, reverse=True):
    return {k: v for k, v in sorted(dict_.items(), key=lambda item: item[1], reverse=reverse)}


def printCodons(seq, start, end, color_seqs=None, 
                highligt_colors = ['on_red', 'on_blue', 'on_magenta'],
                offset=0):
    '''
    Args:
        seq (TYPE): sequence to print in codons.
        start (TYPE): start of the printing.
        end (TYPE): end of the printing.
        color_seqs (TYPE, optional): sequences that will be searched for and 
                                     highlighted with different colors, if found
        highligt_colors: default colors
    '''
    codon_str = ''
    pos_strings = ['' for i in  range(len(str(end)))]
    
    # find matches for color_seqs if specified
    pos_and_color = {}
    pos_and_color['color'] = []
    pos_and_color['idx'] = []
    
    if color_seqs is not None:
        # build index for each provided sequence
        for i, c_seq in enumerate(color_seqs):
            pos_and_color['color'].append(highligt_colors[i])
            
            idx_tmp = findAll(c_seq, seq)
            if len(idx_tmp) == 0:
                continue
            idx_range = []
            for idx_ in idx_tmp:
                idx_range.append(np.arange(idx_,idx_+len(c_seq)))
            
            
            pos_and_color['idx'].append(np.hstack(idx_range))
    
    # for i in range(start,end+1, 3):
    for i in range(start,end+1):
        if color_seqs is None:
            codon_str += seq[i]
        else:
            # for codon_pos in range(3):
            col_tmp = None
            # look for potential coloring
            for indexes, col in zip(pos_and_color['idx'], pos_and_color['color']):
                if i in indexes:
                    col_tmp = col
                    break
            
            if col_tmp is None:
                codon_str += seq[i]
            else:
                codon_str += colored(seq[i] , on_color=col_tmp, attrs=['bold'])
        
        # insert space in between each codons
        if (i-2-offset) % 3 == 0:
            codon_str += ' '
            
        # create string with positions
        for j in range(len(pos_strings)):
            # write out the number under the start of the codon
            if i % 3 == 0:
                if j < len(str(i)):
                    pos_strings[j] += str(i)[j]
                else:
                    pos_strings[j] += ' '                    
            else:
                pos_strings[j] += ' '
            # make an ekstra space if the end of the codon
            if (i-2-offset) % 3 == 0:
                pos_strings[j] += ' '


    # print the results        
    print(codon_str)
    for i in range(len(pos_strings)):
        print(pos_strings[i])
        
def invertSeq(seq):
    inv_dict = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G'
        }
    
    return ''.join([inv_dict[bp] for bp in seq])

file_name = '/media/maltejensen/SSD_EVO_970/phage_data/data'

good_seqs = []
good_qa = []

discarded_seq = {}
'''
Description:
        {'seq': 48 base sequence,
        'qa': quality score of seq,
        'TAG': True if TAG is present,
        'TAG_idx': The index of the TAG, if present,
        'GGGS': True if seq HAS GGGS at the end position,
        'len': length of the whole original sequence,
        'phred_thres': whether below any was below phred threshold 
        }
'''
discard_statistics = {'TAG_bool':0,
                      'GGGS_bool': 0,
                      'phred_thres_bool': 0}

multiple_inserts = []
read = []

phage_start = 127
phage_end = 175 
phred_thres = 30

discard_counter = 0


upstream_seq =   'GGTACC'
downstream_seq = 'CGGCCG'
upstream_start = 105
downstream_start = 173

all_seq = []

t0 = time.time()

with open(file_name, 'r', encoding='ascii') as f:
    for i, line in enumerate(f):
        if i%400000 == 0 and i!=0:
            print('{} reads. {} reads/s'.format(int(i/4), int(100000/(time.time()-t0))) )
            t0 = time.time()
        
        # 0 = header, 1 = sequence, 2 = +, 3 = quality
        read.append(line.strip())
        
        if (i+1) % 4 == 0 and i != 1:    
            
            # all_qa.append(ascii2phred(read[3]))
            
            
            # all_seq.append(findAll(upstream_seq, read[1]))
            # read = []
            # continue
        
            seq_len = len(read[1]) 
            # len_bool = seq_len == 251
            # extract the areas coding for the peptide
            seq_tmp = read[1][phage_start:phage_end]
            qa_tmp = read[3][phage_start:phage_end]
            qa_tmp = ascii2phred(qa_tmp)
            
            if read[1].find(upstream_seq) != upstream_start or read[1].find(downstream_seq) != downstream_start:
                multiple_inserts.append(read[1])
            
            # make checks for exclusion
            TAG_bool, TAG_idx = hasCodon(seq_tmp, 'TAG')
            GGGS_bool = hasGGGSEnd(seq_tmp, GGGS_seqs)
            phred_thres_bool = qa_tmp.mean() < phred_thres
            
            bool_tmp = {
                'TAG_bool': TAG_bool,
                'GGGS_bool': GGGS_bool,
                'phred_thres_bool': False
                }
            
            # check for any errors
            if any(bool_tmp.values()):
                discarded_seq[discard_counter] = {'seq': seq_tmp,
                                                  'qa': qa_tmp,
                                                  'TAG': TAG_bool,
                                                  'TAG_idx': TAG_idx,
                                                  'GGGS': GGGS_bool,
                                                  'len': seq_len,
                                                  'phred_thres': phred_thres_bool                                                  
                    }
                discard_counter += 1
                
                for k, v in bool_tmp.items():
                    if v:
                        discard_statistics[k] += 1
            else:
                good_seqs.append(seq_tmp)
                good_qa.append(qa_tmp)
            
            # empty list for next read
            read = []
#%%
total_n = len(good_seqs) + len(discarded_seq)

dist_1 = 25

print('Report:')
print('{:<{}}{:<10.4e}'.format('Total reads:', dist_1, total_n))
print('{:<{}}{:<15.4e}{:.2f}%'.format('Good reads:', dist_1, len(good_seqs), len(good_seqs)/total_n*100))
print('{:<{}}{:<15.4e}{:.2f}%'.format('Discarded reads:', dist_1, len(discarded_seq), len(discarded_seq)/total_n*100))
print('{:<{}}{:<15.4e}{:.2f}%'.format('Wrongly placed seqs:', dist_1, len(multiple_inserts), len(multiple_inserts)/total_n*100))
for k, v in discard_statistics.items():
    print('{:<{}}{:<15.4e}{:.2f}%'.format(k+':', dist_1, v, v/total_n*100))
#%%

a = 20

print('{:<{}}{}'.format('test name:', a, 200))

#%%
'''
Find location of up and downstream that did not match correctely
'''
upstream_loc = []
downstream_loc = []

upstream_count = {}
downstream_count = {}

upstream_loc_count = {}
downstream_loc_count = {}

upstream_downstream_distance = {}

suspected_frameshifts = []

specific_downstream = []

for m in tqdm(multiple_inserts):
    u_loc = findAll(upstream_seq, m)
    d_loc = findAll(downstream_seq, m)
    
    # append location 
    upstream_loc.append(u_loc)
    downstream_loc.append(d_loc)
    
    # count numbers of copies found in string
    if len(u_loc) not in upstream_count:
        upstream_count[len(u_loc)] = 1
    else: 
        upstream_count[len(u_loc)] += 1
   
    if len(d_loc) not in downstream_count:
        downstream_count[len(d_loc)] = 1
    else: 
        downstream_count[len(d_loc)] += 1
    
    # count placement of the copies 
    for ul in u_loc:
        if ul not in upstream_loc_count:
            upstream_loc_count[ul] = 1
        else: 
            upstream_loc_count[ul] += 1
       
    for ud in d_loc:
        if ud not in downstream_loc_count:
            downstream_loc_count[ud] = 1
        else: 
            downstream_loc_count[ud] += 1
    
    if len(u_loc) == 1 and len(d_loc) == 1:
        # dist = (d_loc[0] - u_loc[0]) - (len(upstream_seq) + 16)
        dist = d_loc[0] - u_loc[0] 
        if dist not in upstream_downstream_distance:
            upstream_downstream_distance[dist] = 1
        else:
            upstream_downstream_distance[dist] += 1
   
        if u_loc[0] == upstream_start and d_loc[0] != downstream_start:
            suspected_frameshifts.append(m)
        
        
        
        GTGA = findAll('GTGA', m)
        
        
        if len(GTGA) > 0:
            # if (u_loc[0]-GTGA[0]) - 32 in [-2,-1,1,2] and (d_loc[0]-GTGA[0])-100 in [-2,-1,1,2] and GTGA[0] in [71, 72, 73, 74, 75]: # 125 for the missing insert
            if (u_loc[0]-GTGA[0]) - 32 in [-2,-1,1,2] and (d_loc[0]-GTGA[0])-100 in [-2,-1,1,2] and GTGA[0] in [71, 72, 73, 74, 75]: # 125 for the missing insert
                specific_downstream.append(m)
                
      
        # if len(GTGA) > 0:
        #     if u_loc[0]-GTGA[0] == 32 and d_loc[0]-GTGA[0] == 100 and GTGA[0] in [71, 72, 73, 74, 75]: # 125 for the missing insert
        #         specific_downstream.append(m)
                
      
        
      
        # GTG_bool = any([num in findAll('GTG', m) for num in [71, 72, 73, 74, 75]])
        
        
        
        # GTGA = findAll('GTGA', m)
        
        
        # if len(GTGA) > 0:
        #     GTGA_bool = GTGA[0] == 73 
        #     u_loc_bool = u_loc in [103,104,106,107]
        #     d_loc_bool = d_loc in [171,172,174,175]
            
        #     if all([GTGA_bool, u_loc_bool, d_loc_bool]):
        #         specific_downstream.append(m)
            
        
        # GTG_bool = any([num in findAll('GTG', m) for num in [71, 72, 74, 75]])
        # u_loc_bool = u_loc in [103,104,106,107]
        # d_loc_bool = d_loc in [171,172,174,175]
        
        # if any([GTG_bool, u_loc_bool, d_loc_bool]):
        #     specific_downstream.append(m)
        
        
            
            
upstream_count = dict(sorted(upstream_count.items()))
downstream_count = dict(sorted(downstream_count.items()))
upstream_loc_count = dict(sorted(upstream_loc_count.items()))
downstream_loc_count = dict(sorted(downstream_loc_count.items()))
upstream_downstream_distance = dict(sorted(upstream_downstream_distance.items()))

'''
Many downstream on position 125: 'TCT' just before the insert is missing, and the 
downstream part is directly after

For downstream at 170, a whole codon is missing from the insert


'''
#%%
for i in range(20):
    # printCodons(specific_downstream[i], 69,178, color_seqs=[upstream_seq, downstream_seq, 'TTCTATTCTCACTCT'])
    printCodons(specific_downstream[i], 69,178, color_seqs=[upstream_seq, downstream_seq, 'GTGA'], offset=1)

    
#%%
print('Number of upstream seqs')
for k, v in upstream_count.items():
    print('{:<5}{:<10}{:.2f}%'.format(str(k)+':',v, v/np.sum(list(upstream_count.values()))*100 ))

print('\nNumber of downstream seqs')
for k, v in downstream_count.items():
    print('{:<5}{:<10}{:.2f}%'.format(str(k)+':',v, v/np.sum(list(upstream_count.values()))*100 ))
#%%
print('Placements of upstream seqs')
for k, v in upstream_loc_count.items():
    print('{:<5}{:<12}{:.3f}%'.format(str(k)+':',v, v/np.sum(list(upstream_loc_count.values()))*100 ))
#%%
print('Placements of downstream seqs')
for k, v in downstream_loc_count.items():
    print('{:<5}{:<12}{:.3f}{:<9}{:.3f}%'.format(str(k)+':',v, 
                                                v/np.sum(list(downstream_loc_count.values()))*100, 
                                                '%', v/total_n*100 ))
#%%
print('Distances between upstream and downstream seqs')
for k, v in upstream_downstream_distance.items():
    print('{:<5}{:<12}{:.3f}%'.format(str(k)+':',v, v/np.sum(list(upstream_downstream_distance.values()))*100 ))

#%%
stacked_qa = np.hstack(all_qa)
max_qa = stacked_qa.max()
#%%
plt.figure()
plt.hist(stacked_qa, bins = np.arange(max_qa+2), density=True)
plt.ylabel('Fraction')
plt.xlabel('Phred Score')
plt.show()
#%%

v_stacked_qa = np.vstack(all_qa)
#%%


plt.figure(figsize=(15,10))
# plt.errorbar(np.arange(1,v_stacked_qa.shape[1]+1), v_stacked_qa[:1000,:].mean(axis=0), 
#               yerr=v_stacked_qa[:1000,:].std(axis=0), fmt = '.', barsabove=False, capsize=3, ecolor='gray')
plt.errorbar(np.arange(1,v_stacked_qa.shape[1]+1), v_stacked_qa.mean(axis=0), 
              yerr=v_stacked_qa.std(axis=0), fmt = '.', barsabove=False, capsize=3, ecolor='gray')
# plt.axvline(x=36.5)
# plt.axvline(x=48.5)
plt.axvline(x=105.5)
plt.axvline(x=180.5)
plt.axvline(x=phage_start, color='green')
plt.axvline(x=phage_end-3, color='green')
plt.grid()
plt.tight_layout()
plt.show()    
#%%
plt.figure()
plt.boxplot(v_stacked_qa)
plt.show()
#%%
plt.xlabel('Nucleotide Position')
plt.ylabel('Phred Score')

#%%

tag_pos = []

for k, v in discarded_seq.items():
    if v['TAG']:
        tag_pos.append(v['TAG_idx'])


print('number of TAG: {:.5e}'.format(len(tag_pos)))
#%%
plt.figure()
plt.hist(np.hstack(tag_pos), bins=np.arange(1,18), density=True)

#%%
pos_len = []

for pos in tag_pos:
    pos_len.append(len(pos))

pos_len = np.array(pos_len)

for len_ in range(1,pos_len.max()):
    
    num_len = (pos_len == len_).sum()
    print('{} TAG in sequence:{:>15}{:>10.2f}%'.format(len_, num_len, num_len/len(pos_len)*100 ))

#%%
'''
Dont have GGGS
'''
gggs_pos = []

for k, v in discarded_seq.items():
    if not v['GGGS']:
        gggs_pos.append(v['seq'])

print('number of GGGS: {:.5e}'.format(len(gggs_pos)))
#%%
'''
Both TAG and dont have GGGS
'''
gggs_and_TAG = []

for k, v in discarded_seq.items():
    if not v['GGGS'] and v['TAG']:
        gggs_and_TAG.append(v)

print('number of TAG and GGGS: {:.5e}'.format(len(gggs_and_TAG)))

#%%
discarded_qa = []

for k, v in discarded_seq.items():
    discarded_qa.append(v['qa'])

stacked_discarded_qa = np.hstack(discarded_qa)

max_qa = stacked_discarded_qa.max()

plt.figure()
plt.hist(stacked_discarded_qa, bins = np.arange(max_qa+1), density=True, label='Discarded')
plt.hist(stacked_qa, bins = np.arange(max_qa+1), density=True,alpha=0.7, label='Good')
plt.ylabel('Fraction')
plt.xlabel('Phred Score')
plt.legend()
plt.show()
#%%

v_stacked_discarded_qa = np.vstack(discarded_qa)
#%%

plt.errorbar(np.arange(1,v_stacked_discarded_qa.shape[1]+1), v_stacked_discarded_qa.mean(axis=0), 
              yerr=v_stacked_discarded_qa.std(axis=0), fmt = '.', barsabove=False, capsize=3, ecolor='gray')
plt.axvline(x=36.5)
plt.axvline(x=48.5)
plt.grid()
plt.show()    

#%%

all_qa = np.hstack((stacked_qa, stacked_discarded_qa))

print('Expected P of mistake: {:.2e}'.format( (10.0**(-all_qa/10)).mean()) )

# Expected P of mistake: 5.02e-03
#%%
from scipy.stats import binom

k = 1
n = 12

print('{} copies: {:.2e}'.format(k, binom.pmf(k, 12, 1.04e-4)))
#%%
'''
Inspect the distribution of the last 4 wrongs codons (not GGGS)
'''
# most prevalent: GGT GGA GGT TCG -> GGGS (511,138)
# GGA GGA GGT TCG

only_end = {}

for i, seq_ in enumerate(gggs_pos):
    if seq_[-12:] not in only_end:
        only_end[seq_[-12:]] = 1
    else:
        only_end[seq_[-12:]] += 1


only_end = {k: v for k, v in sorted(only_end.items(), key=lambda item: item[1], reverse=True)}

for k, v in only_end.items():
    print('{}:{:<15}'.format(k,v))
#%%
'''
Inspect the distribution of the last 4 wrongs codons (not GGGS)
'''
# most prevalent: GGT GGA GGT TCG -> GGGS (511,138)
# GGT GGA GGT TCG
only_end_good = {}

for i, seq_ in enumerate(good_seqs):
    if seq_[-12:] not in only_end_good:
        only_end_good[seq_[-12:]] = 1
    else:
        only_end_good[seq_[-12:]] += 1


only_end_good = {k: v for k, v in sorted(only_end_good.items(), key=lambda item: item[1], reverse=True)}

for k, v in only_end_good.items():
    print('{}:{:<15}'.format(k,v))

#%%
'''
Save the data
'''
with open('/media/maltejensen/SSD_EVO_970/phage_data/test_reads.txt', 'r', encoding='ascii') as f:
    line = f.readline()
    
    
    while line:
        print(line.strip())
        line = f.readline()
    
#%%
'''
Preparation for weblogo
'''
# output name to use for weblogo
save_name = '/media/maltejensen/SSD_EVO_970/phage_data/extracted_sequences_2.txt'

with open(save_name, 'w') as f:
    for i, seq in enumerate(good_seqs):
        f.write(seq[:-12]+'\n')


#%%
import weblogo

'''
Create weblogo
'''

with open(save_name, 'r') as f:    
    # seqs = weblogo.seq_io.read_seq(f)
    # seqs = weblogo.seq_io.array_io.read(f, alphabet=weblogo.seq.Alphabet('ATGC'))
    seqs = weblogo.logo.read_seq_data(f)


logodata = weblogo.LogoData.from_seqs(seqs)


logooptions = weblogo.LogoOptions()
logooptions.title = "A Logo Title"
# logooptions.resolution = 300
logooptions.xaxis_label = 'Position'
logooptions.show_fineprint = False
logooptions.unit_name = 'probability'
logoformat = weblogo.LogoFormat(logodata, logooptions)
# logoformat.color_scheme = ('red', 'green')
logoformat.color_scheme = weblogo.colorscheme.nucleotide

# eps = weblogo.eps_formatter(logodata, logoformat)
# png = weblogo.png_formatter(logodata, logoformat)
png = weblogo.png_print_formatter(logodata, logoformat)

# logodata = weblogo.LogoData.from_seqs(seqs)


import matplotlib.pyplot as plt
import io 
import matplotlib.image as mpimg

fp = io.BytesIO(png)

with fp:
    img = mpimg.imread(fp, format='jpeg')


fig, axx = plt.subplots(1, figsize=(15,5))
axx.imshow(img)
# axx.xaxis.set_ticklabels([])
# axx.yaxis.set_ticklabels([])
# axx.set_xticks([])
# axx.set_yticks([])
axx.set_axis_off()


# plt.imshow(img)
plt.show()
#%%
'''
Make heatmap for nucleotides 
'''
from tqdm import tqdm

def OneHotEncodeDNA(seq, map_):


    num_col = len(map_)
    one_hot = np.zeros(num_col*len(seq), dtype=np.int8)
    
    for i in range(len(seq)):
        one_hot[i*num_col + map_[seq[i]]] = 1

    return one_hot

dna_map_ = {
    'A': 0,
    'T': 1,
    'G': 2,
    'C': 3
    }

one_hot_dna = []

N_counter = 0

for seq in tqdm(good_seqs):
    if 'N' in seq:
        N_counter += 1
        continue
    
    one_hot_dna.append(OneHotEncodeDNA(seq[:36], dna_map_))
    
one_hot_dna = np.vstack(one_hot_dna)

#%%

summed_nucleotides = one_hot_dna.sum(axis=0)/one_hot_dna.shape[0]
nucleotides =list( dna_map_.keys())

for i in range(one_hot_dna.shape[1]//4):
    for j in range(4):
        print('{}: {:.5f}'.format('('+str(i+1)+')'+nucleotides[j], summed_nucleotides[i*4+j]))
    print('')

#%%
'''
Heatmap if nucleotides
'''
heatmap = one_hot_dna.sum(axis=0)/one_hot_dna.shape[0]
heatmap = heatmap.reshape(4,-1, order='F')


heatmap = heatmap-np.tile([[0.25,0.25, 0], [0.25,0.25, 0.5], [0.25,0.25, 0.5], [0.25,0.25, 0]],12)
heatmap *= 100

#%%
from mpl_toolkits.axes_grid1 import make_axes_locatable
abs_max = np.max(np.abs(heatmap))


fig, axx = plt.subplots(1, figsize=(13,3))

img = axx.imshow(heatmap, cmap='bwr', vmin=-abs_max, vmax=abs_max)
axx.set_yticks(np.arange(len(dna_map_)))
axx.set_yticklabels(list(dna_map_.keys()))
axx.set_xlabel('Position')


# axx.set_xticks(np.arange(0,36,3))
# axx.set_xticklabels(list(np.arange(0,36,3)+1 ))
axx.set_xticks(np.arange(36))
axx.set_xticklabels(np.arange(36)+1)

divider = make_axes_locatable(axx)
cax = divider.append_axes("bottom", size="15%", pad=0.5)
cbar = plt.colorbar(img, cax=cax, label='Difference from expected [%]', orientation ='horizontal')


plt.show()
#%%

aa_sub = df['AA'][pred_ == 3]

df_dict = {i+1: [] for i in range(num_seq)}

for aa in aa_sub:
    for i, a in enumerate(list(aa)[:-4], start=1):
        df_dict[i].append(a)

df_dict = pd.DataFrame(df_dict)
        
#%%
from mlxtend.frequent_patterns import apriori, association_rules

# cluster_subset = one_hot[pred_ == 2-1, :]
cluster_subset = one_hot

sub_df = pd.DataFrame(cluster_subset, dtype=bool)



#%%
import time

t0 = time.time()

frq_items = apriori(sub_df, min_support = 0.02, use_colnames = True, verbose=2)
frq_items = frq_items.sort_values(by='support', ascending=False)
frq_items = frq_items.reset_index(drop=True)

print('Time: {:.1f} min'.format((time.time()-t0)/60 ))

#%%

aa_freq = {k: v/sum(num_codons.values()) for k, v in num_codons.items()}

perc = sub_df.shape[0]/one_hot.shape[0]

def getAAPos(one_hot, seq_map_rev):
    aa = seq_map_rev[one_hot%len(seq_map_rev)]
    pos = one_hot//len(seq_map_rev)
    return aa, pos

item_str_list = []
support_list = []
adj_support = []
expected_list = []
diff_list = []
diff_frac_list = []
num_aa_list = []

# print('{:<20}{:<15}{:<15}{:<15}{}'.format('Rule', 'Support', 'Adj. support', 
                                    # 'Expected', 'Diff'))
for i, row in frq_items.iterrows():
    items = list(row['itemsets'])
    items.sort()
    item_str = ''
    
    if len(items)>0:
        prob = 1
        
        for it in items:
            aa_ = seq_map_rev[it%num_groups]
            prob *= aa_freq[aa_]
            item_str += '{}({}) '.format(aa_, it//num_groups+1)
    
        
        item_str_list.append(item_str)
        support_list.append(row['support']*100)
        adj_support.append(row['support']*perc*100)
        expected_list.append(100*prob)
        diff_list.append((row['support']*perc - prob)*100)
        diff_frac_list.append((row['support']*perc / prob))
        num_aa_list.append(len(items))

    
    #     print('{:<20}{:<15.3f}{:<15.3f}{:<15.3f}{:.3f}'.format(item_str, 
    #                 row['support']*100, row['support']*perc*100, 
    #                 100*prob, (row['support']*perc - prob)*100))
    # # if i == 30:
    #     break

frq_items_proc = pd.DataFrame({
    'Rule': item_str_list,
    'Support': support_list, 
    'Adj. support': adj_support, 
    'Expected': expected_list, 
    'Diff': diff_list,
    'Diff frac.': diff_frac_list,
    'Num AA': num_aa_list
    })

frq_items_proc = frq_items_proc.sort_values(by=['Num AA','Support'], ascending=False)
# frq_items_proc = frq_items_proc.sort_values(by=['Num AA','Diff frac.'], ascending=False)
frq_items_proc = frq_items_proc.reset_index(drop=True)
#%%

# def mySort(a):
#     print(a)
    # return len(a)

# def mySort(a,b):
#     return len(a['antecedents']) > len(b['antecedents'])

# # Collecting the inferred rules in a dataframe
rules = association_rules(frq_items, metric ='confidence', min_threshold = 0.8)


for i, row in rules.iterrows():
    ante = list(row['antecedents'])
    conse = list(row['consequents'])
    ante.sort()
    conse.sort()
    ante_str = ''
    conse_str = ''
    
    for an in ante:
        aa_ = seq_map_rev[an%num_groups]
        ante_str += '{}({}) '.format(aa_, an//num_groups+1)
    
    rules['antecedents'][i] = ante_str
   
    for co in conse:
        aa_ = seq_map_rev[co%num_groups]
        conse_str += '{}({}) '.format(aa_, co//num_groups+1)
    
    rules['consequents'][i] = conse_str
   

# rules = rules.sort_values(by='antecedents', key=mySort, axis=0)

# rules = rules.sort_values(['confidence', 'lift '], ascending =[False, False])
# print(rules.head())
