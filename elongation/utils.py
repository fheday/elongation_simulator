'''
Utility functions
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
from collections import Counter

def get_number_elongating_ribosomes_per_time(sim):
    return [len(entry) for entry in sim.ribosome_positions_history], np.cumsum(sim.dt_history)


'''
Returns a tuple where the first element is a vector with the
enlogation duration of the ribosomes that terminated in the simulation, and
the second element is a vector with the iteration where such ribosomes
started enlogating.
'''
def get_codon_average_occupancy(sim):
    # sim.updateRibosomeHistory()
    return sim.getEnlogationDuration()


canonical_bases = ['A', 'C', 'U', 'G']
bases = ['&', 'N', '3', '1', 'P', '~', 'I', '#', '?'] + canonical_bases


def first_WC(codon, tRNA):
    codon_nt = codon[0] # first character
    tRNA_nt = tRNA[2] # thrird character
    return (codon_nt == 'A' and tRNA_nt == 'U') or (codon_nt == 'C' and tRNA_nt == 'G') or \
           (codon_nt == 'G' and tRNA_nt == 'C') or (codon_nt == 'U' and tRNA_nt == 'A')

def second_WC(codon, tRNA):
    codon_nt = codon[1] # second character
    tRNA_nt = tRNA[1] # second character
    return (codon_nt == 'A' and tRNA_nt == 'U') or (codon_nt == 'C' and tRNA_nt == 'G') or \
           (codon_nt == 'G' and tRNA_nt == 'C') or (codon_nt == 'U' and tRNA_nt == 'A')

def third_WC(codon, tRNA, gene):
    codon_nt = codon[2] # thrid character
    tRNA_nt = tRNA[0] # thrird character
    indexes = []
    for i in range(4):
        indexes.append(np.argwhere(np.array(gene[i*13:(i+1) * 13])==1)[:,0])
    return (codon_nt == 'A' and tRNA_nt in ['U'] + [bases[i] for i in indexes[0]]) or (codon_nt == 'C' and tRNA_nt in ['G'] + [bases[i] for i in indexes[1]]) or \
           (codon_nt == 'G' and tRNA_nt in ['C'] + [bases[i] for i in indexes[2]]) or (codon_nt == 'U' and tRNA_nt in ['A'] + [bases[i] for i in indexes[3]])

def first_wobble(codon, tRNA, gene):
    codon_nt = codon[0] # first character
    tRNA_nt = tRNA[2] # thrird character
    indexes = []
    for i in range(4):
        indexes.append(np.argwhere(np.array(gene[52+i*4:52 + (i+1) * 4])==1)[:,0])
    return (codon_nt == "A" and tRNA_nt in [canonical_bases[i] for i in indexes[0]]) or (codon_nt == "C" and tRNA_nt in [canonical_bases[i] for i in indexes[1]]) or \
           (codon_nt == "G" and tRNA_nt in [canonical_bases[i] for i in indexes[2]]) or (codon_nt == "U" and tRNA_nt in [canonical_bases[i] for i in indexes[3]])

def third_wobble(codon, tRNA, gene):
    codon_nt = codon[2] # thrid character
    tRNA_nt = tRNA[0] # thrird character
    indexes = []
    for i in range(4):
        indexes.append(np.argwhere(np.array(gene[68 + i * 13:68 + (i+1) * 13])==1)[:,0])
    return (codon_nt == "A" and tRNA_nt in [bases[i] for i in indexes[0]]) or (codon_nt == "C" and tRNA_nt in [bases[i] for i in indexes[1]]) or \
           (codon_nt == "G" and tRNA_nt in [bases[i] for i in indexes[2]]) or (codon_nt == "U" and tRNA_nt in [bases[i] for i in indexes[3]])

# make.matrix
def make_matrix(tRNAs, codons, gene, verbose=False):
    cognate_WC_matrix = np.zeros((len(tRNAs.anticodon), len(codons.codon)))
    cognate_wobble_matrix = np.zeros((len(tRNAs.anticodon), len(codons.codon)))
    nearcognate_matrix = np.zeros((len(tRNAs.anticodon), len(codons.codon)))
    
    if verbose:
        print("Populating WC matrix...")
    # populate cognate WC matrix if WC criteria matched
    for n in range(len(tRNAs.anticodon)):
        for m in range(len(codons.codon)):
            codon = codons.codon[m]
            anticodon = tRNAs.anticodon[n]
            if second_WC(codon, anticodon) and first_WC(codon, anticodon) and third_WC(codon, anticodon, gene):
                cognate_WC_matrix[n, m] = 1
    if verbose:
        print("done.")
        print("Populating wobble matrix...")
        
    
    #populate cognate wobble matrix if wobble criteria matched, amino acid is correct, and WC matrix entry is 0
    #if incorrect amino acid, assign to near-cognates
    for n in range(len(tRNAs.anticodon)):
        for m in range(len(codons.codon)):
            if cognate_WC_matrix[n,m] == 0 and second_WC(codons.codon[m],tRNAs.anticodon[n]) and\
            first_WC(codons.codon[m],tRNAs.anticodon[n]) and third_wobble(codons.codon[m],tRNAs.anticodon[n], gene):
                if tRNAs["three.letter"][n] == codons["three.letter"][m]:
                    cognate_wobble_matrix[n,m] = 1
                else:
                    nearcognate_matrix[n,m] = 1  

    if verbose:
        print('done.')
        print('Populating nearcognate matrix...')

    #populate near-cognate matrix if wobble criteria matched and wobble and WC matrix entries are 0
    for n in range(len(tRNAs.anticodon)):
        for m in range(len(codons.codon)):
            if (cognate_WC_matrix[n,m] == 0 and cognate_wobble_matrix[n,m] == 0 and second_WC(codons.codon[m],tRNAs.anticodon[n]) and \
                first_wobble(codons.codon[m],tRNAs.anticodon[n], gene) and third_wobble(codons.codon[m],tRNAs.anticodon[n], gene)):
                nearcognate_matrix[n,m] = 1

    if verbose:
        print('done.')

    return {"cognate.wc.matrix":cognate_WC_matrix, "cognate.wobble.matrix":cognate_wobble_matrix, "nearcognate.matrix":nearcognate_matrix}

def is_valid_matrix(matrices_dict, verbose=False):
    #Sanity checks

    #Check whether any tRNA:codon combination is assigned 1 in more than one table (this should not occur)
    if len(matrices_dict.keys()) == 0:
        if verbose: print("no keys in matrices_dict.")
        return False
    testsum = np.zeros((matrices_dict[list(matrices_dict.keys())[0]].shape))
    for k in matrices_dict.keys():
        testsum += matrices_dict[k]
#     testsum = cognate_WC_matrix + cognate_wobble_matrix + nearcognate_matrix
    if np.any(testsum>1):
        if verbose: print('Warning: multiple relationships for identical tRNA:codon pairs detected.')
        return False
    elif np.any(np.sum(testsum, axis=0) == 0):
        if verbose: print('Warning: codon without decoding tRNA detected.')
        return False
    elif verbose:
        print('No threesome errors detected.')
    return True

def plot_matrix(codons, tRNAs, matrices_dict):
    colours=['g', 'y', 'r']
    labels = list(matrices_dict.keys())
    i = 0
    plt.figure(figsize=(25,15))
    plt.grid(True)
    for k in matrices_dict.keys():
        c = np.argwhere(matrices_dict[k] == 1)#[:,1]
        plt.plot(c[:,1], c[:,0], colours[i] + 's', label=labels[i])
        i +=1
    plt.xticks(range(len(codons.codon)), codons.codon, rotation = 45)
    plt.yticks(range(len(tRNAs.anticodon)), tRNAs.anticodon)
    plt.legend()
    
def make_concentrations(tRNAs, codons, matrices, total_tRNA=190):

    WCcognate = matrices["cognate.wc.matrix"]
    wobblecognate = matrices["cognate.wobble.matrix"]
    nearcognate = matrices["nearcognate.matrix"]
    if "gene.copy.number" not in tRNAs.columns:
        print("tRNA list must contain columns with abundance information headed either 'gene.copy.number' or 'experimental.abundance'.")
        return pd.DataFrame()
    
    # construct empty results dataframe
    tRNA_concentrations = pd.DataFrame(codons[codons.columns[:2]])
    tRNA_concentrations["WCcognate.conc"] = 0.0
    tRNA_concentrations["wobblecognate.conc"] = 0.0
    tRNA_concentrations["nearcognate.conc"] = 0.0
    
    #calculate a conversion factor to convert the abundance factor to a molar concentration
    conversion_factor = total_tRNA / tRNAs["gene.copy.number"].sum() * 1e-6
    
    #go through the WCcognates matrix and for each entry of 1 add the abundance of the tRNA from the abundance table to the concentration table
    for n in range(len(codons.codon)):
        for m in range(len(tRNAs.anticodon)):
            if WCcognate[m, n] == 1:
                tRNA_concentrations.loc[n, "WCcognate.conc"] = tRNA_concentrations["WCcognate.conc"][n] + (tRNAs["gene.copy.number"][m]*conversion_factor)
    
    #ditto for wobblecognate
    for n in range(len(codons.codon)):
        for m in range(len(tRNAs.anticodon)):
            if wobblecognate[m, n] == 1:
                tRNA_concentrations.loc[n, "wobblecognate.conc"] = tRNA_concentrations["wobblecognate.conc"][n] + (tRNAs["gene.copy.number"][m]*conversion_factor)

    #ditto for nearcognates
    for n in range(len(codons.codon)):
        for m in range(len(tRNAs.anticodon)):
            if nearcognate[m, n] == 1:
                tRNA_concentrations.loc[n, "nearcognate.conc"] = tRNA_concentrations["nearcognate.conc"][n] + (tRNAs["gene.copy.number"][m]*conversion_factor)
    return tRNA_concentrations
