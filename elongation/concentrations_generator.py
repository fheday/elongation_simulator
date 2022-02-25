#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 02 11:32:41 2019

Re-implementation of the R-code that calculates concentrations
This script is intented to generate new concentrations (WC, wobble, near) with different interpretation
of wobble and near cognates. it is intended to be use in GA algorithm.

@author: heday
"""

from cmath import nan
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

# make.matrix
"""
Given a DataFrame with tRNA concentrations, and another DataFrame with codons information, generates a decoding matrix

TRNAs: DataFrame with the tRNAs: anticodon, gene.copy.number
Codons: DataFrame with the codons to be used
"""
def make_matrix(tRNAs, codons, verbose=False):
    # check if tRNAs have 'anticodon' column
    if 'anticodon' not in tRNAs.columns:
        print('tRNA list must contain a column named "anticodon".')
        return
    # check if codons have 'codon' column
    if 'codon' not in codons.columns:
        print('Codon list must contain a column named "codon".')
        return

    def first_WC(codon, tRNA):
        codon_nt = codon[0] # first character
        tRNA_nt = tRNA[2] # thrird character
        return (codon_nt == "A" and tRNA_nt == "U") or (codon_nt == "C" and tRNA_nt == "G") or \
               (codon_nt == "G" and tRNA_nt == "C") or (codon_nt == "U" and tRNA_nt == "A")

    def second_WC(codon, tRNA):
        codon_nt = codon[1] # second character
        # print(tRNA)
        if tRNA is np.NaN: 
            return False
        tRNA_nt = tRNA[1] # second character
        return (codon_nt == "A" and tRNA_nt == "U") or (codon_nt == "C" and tRNA_nt == "G") or \
               (codon_nt == "G" and tRNA_nt == "C") or (codon_nt == "U" and tRNA_nt == "A")

    def third_WC(codon, tRNA):
        codon_nt = codon[2] # thrid character
        tRNA_nt = tRNA[0] # thrird character
        return (codon_nt == "A" and tRNA_nt in ['U','&','N','3','1','P','~','V','}','S',')','{']) or (codon_nt == "C" and tRNA_nt in ['G','#', 'Q']) or \
               (codon_nt == "G" and tRNA_nt in ['C','B','?', 'M']) or (codon_nt == "U" and tRNA_nt in ['A','I'])

    def first_wobble(codon, tRNA):
        codon_nt = codon[0] # first character
        tRNA_nt = tRNA[2] # thrird character
        return (codon_nt == "A" and tRNA_nt in ["A", "U"]) or (codon_nt == "C" and tRNA_nt in ["G", "A", "U"]) or \
               (codon_nt == "G" and tRNA_nt in ["A", "C", "U"]) or (codon_nt == "U" and tRNA_nt in ["A", "G", "U"])

    def third_wobble(codon, tRNA):
        codon_nt = codon[2] # thrid character
        tRNA_nt = tRNA[0] # thrird character
        return (codon_nt == "A" and tRNA_nt in ['A','U','&','N','3','1','P','~','I']) or (codon_nt == "C" and tRNA_nt in ['G','#','A','U','I']) or \
               (codon_nt == "G" and tRNA_nt in ['C','B','?','A','U','&','N','3','1','P','~', 'V', 'S', ')', '{']) or (codon_nt == "U" and tRNA_nt in ['A','G','U','I','#', 'Q', 'V'])

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
            if second_WC(codon, anticodon) and first_WC(codon, anticodon) and third_WC(codon, anticodon):
                cognate_WC_matrix[n, m] = 1
    if verbose:
        print("done.")
        print("Populating wobble matrix...")
        
    
    #populate cognate wobble matrix if wobble criteria matched, amino acid is correct, and WC matrix entry is 0
    #if incorrect amino acid, assign to near-cognates

    for n in range(len(tRNAs.anticodon)):
        for m in range(len(codons.codon)):
            if cognate_WC_matrix[n,m] == 0 and second_WC(codons.codon[m],tRNAs.anticodon[n]) and\
            first_WC(codons.codon[m],tRNAs.anticodon[n]) and third_wobble(codons.codon[m],tRNAs.anticodon[n]):
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
                first_wobble(codons.codon[m],tRNAs.anticodon[n]) and third_wobble(codons.codon[m],tRNAs.anticodon[n])):
                nearcognate_matrix[n,m] = 1

    if verbose:
        print('done.')

    #Sanity checks

    #Check whether any tRNA:codon combination is assigned 1 in more than one table (this should not occur)

    testsum = cognate_WC_matrix + cognate_wobble_matrix + nearcognate_matrix
    if np.any(testsum>1):
      print('Warning: multiple relationships for identical tRNA:codon pairs detected.')
      return {}
    elif verbose:
      print('No threesome errors detected.')

    return {"cognate.wc.matrix":cognate_WC_matrix, "cognate.wobble.matrix":cognate_wobble_matrix, "nearcognate.matrix":nearcognate_matrix}

"""
Plots the pairing matrices.
"""
def plot_matrix(matrices_dict, tRNAs, codons):
    colours=['g', 'y', 'r']
    labels = list(matrices_dict.keys())
    i = 0
    plt.figure(figsize=(25,15))
    plt.grid(True)
    xmax = 0
    for k in matrices_dict.keys():
        c = np.argwhere(matrices_dict[k] == 1)#[:,1]
        # display(c)
        plt.plot(c[:,1], c[:,0], colours[i] + 's', label=labels[i])
        i +=1
    plt.xticks(range(len(codons.codon)), codons.codon, rotation = 45)
    plt.yticks(range(len(tRNAs.anticodon)), tRNAs.anticodon)
    plt.legend()
    
"""
Given a tRNA matrix, and the decoding matrix, generates a concentrations DataFrame.

TRNAs: DataFrame with the tRNAs: anticodon, gene.copy.number
Matrices: pairing matrices generated by make_matrix
Codons: DataFrame with the codons to be used
use_gene_copy_number: if True, will use tRNAs column ‘gene.copy.number’, otherwise will use ‘experimental.abundance’
total_Trna: default value is 190 (micro moles).
"""
def make_concentrations(tRNAs, matrices, codons, use_gene_copy_number=True, total_tRNA=190):
    #check if directory exists
    # if (not os.path.isdir(directory)):
    #     print('Directory ', directory, 'not found.')

    # #check if codons.csv exits
    # codons = None
    # if os.path.isfile(os.path.join(directory, "codons.csv")):
    #     codons = pd.read_csv(os.path.join(directory, "codons.csv"))
    WCcognate = matrices["cognate.wc.matrix"]
    wobblecognate = matrices["cognate.wobble.matrix"]
    nearcognate = matrices["nearcognate.matrix"]
    # if "gene.copy.number" not in tRNAs.columns:
    #     print("tRNA list must contain columns with abundance information headed either 'gene.copy.number' or 'experimental.abundance'.")
    #     return pd.DataFrame()
    
    concentration_col_name = 'gene.copy.number'
    if "gene.copy.number" not in tRNAs.columns and "experimental.abundance" not in tRNAs.columns:
        print("tRNA list must contain columns with abundance information headed either 'gene.copy.number' or 'experimental.abundance'.")
        return pd.DataFrame()
    # if "gene.copy.number" in tRNAs.columns and "experimental.abundance" in tRNAs.columns:
    if use_gene_copy_number:
        concentration_col_name = 'gene.copy.number'
    else:
        concentration_col_name = 'experimental.abundance'

    # construct empty results dataframe
    tRNA_concentrations = pd.DataFrame(codons[[codons.columns[0],codons.columns[2]]])
    tRNA_concentrations["WCcognate.conc"] = 0.0
    tRNA_concentrations["wobblecognate.conc"] = 0.0
    tRNA_concentrations["nearcognate.conc"] = 0.0
    
    #calculate a conversion factor to convert the abundance factor to a molar concentration
    print('using: '+concentration_col_name)
    conversion_factor = np.float64(total_tRNA / np.float64(tRNAs[concentration_col_name].sum()) * 1e-6)
    
    #go through the WCcognates matrix and for each entry of 1 add the abundance of the tRNA from the abundance table to the concentration table
    for n in range(len(codons.codon)):
        for m in range(len(tRNAs.anticodon)):
            if WCcognate[m, n] == 1:
                tRNA_concentrations.loc[n, "WCcognate.conc"] = tRNA_concentrations["WCcognate.conc"][n] + (tRNAs[concentration_col_name][m]*conversion_factor)
    
    #ditto for wobblecognate
    for n in range(len(codons.codon)):
        for m in range(len(tRNAs.anticodon)):
            if wobblecognate[m, n] == 1:
                tRNA_concentrations.loc[n, "wobblecognate.conc"] = tRNA_concentrations["wobblecognate.conc"][n] + (tRNAs[concentration_col_name][m]*conversion_factor)

    #ditto for nearcognates
    for n in range(len(codons.codon)):
        for m in range(len(tRNAs.anticodon)):
            if nearcognate[m, n] == 1:
                tRNA_concentrations.loc[n, "nearcognate.conc"] = tRNA_concentrations["nearcognate.conc"][n] + (tRNAs[concentration_col_name][m]*conversion_factor)
    return tRNA_concentrations

##example of how to use:
# matrices_dict = conc_generator.make_matrix(tRNAs)
# #here we can change the column gene.copy.number in tRNAs
# x=conc_generator.make_concentrations(tRNAs, matrices_dict)

# tRNAs = pd.read_csv('/home/heday/Projects/R3/Native Spike and B117 Kent/HEK293_processed.csv')
# codons = pd.read_csv('/home/heday/Projects/R3/Native Spike and B117 Kent/codons.csv')
# matrix = make_matrix(tRNAs, codons, verbose=True)