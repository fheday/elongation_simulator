#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 02 11:32:41 2019

Re-implementation of the R-code that calculates concentrations
This script is intented to generate new concentrations (WC, wobble, near) with different
interpretation of wobble and near cognates. it is intended to be use in GA algorithm.

@author: heday
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def make_matrix(t_rnas, codons, verbose=False):
    """
    Given a DataFrame with tRNA concentrations, and another DataFrame with codons information,
    generates a decoding matrix

    TRNAs: DataFrame with the tRNAs: anticodon, gene.copy.number
    Codons: DataFrame with the codons to be used
    """
    # check if tRNAs have 'anticodon' column
    if 'anticodon' not in t_rnas.columns:
        print('tRNA list must contain a column named "anticodon".')
        return
    # check if codons have 'codon' column
    if 'codon' not in codons.columns:
        print('Codon list must contain a column named "codon".')
        return

    def first_wc(codon, trna):
        codon_nt = codon[0] # first character
        trna_nt = trna[2] # third character
        return (codon_nt == "A" and trna_nt == "U") or (codon_nt == "C" and trna_nt == "G") or \
               (codon_nt == "G" and trna_nt == "C") or (codon_nt == "U" and trna_nt == "A")

    def second_wc(codon, trna):
        codon_nt = codon[1] # second character
        trna_nt = trna[1] # second character
        return (codon_nt == "A" and trna_nt == "U") or (codon_nt == "C" and trna_nt == "G") or \
               (codon_nt == "G" and trna_nt == "C") or (codon_nt == "U" and trna_nt == "A")

    def third_wc(codon, trna):
        codon_nt = codon[2] # third character
        trna_nt = trna[0] # first character
        return (codon_nt == "A" and trna_nt in ['U','&','3','1','~','N','S',')','{','V','}','P']) or \
                (codon_nt == "C" and trna_nt in ['G','#','W']) or \
                (codon_nt == "G" and trna_nt in ['C','B']) or (codon_nt == "U" and trna_nt in ['A'])

    def first_wobble(codon, trna):
        codon_nt = codon[0] # first character
        trna_nt = trna[2] # third character
        return (codon_nt == "A" and trna_nt in ['A']) or \
                (codon_nt == "C" and trna_nt in ['A','U']) or \
                (codon_nt == "G" and trna_nt in ['A','U']) or \
                (codon_nt == "U" and trna_nt in ['G','U'])

    def third_wobble(codon, trna):
        codon_nt = codon[2] # third character
        trna_nt = trna[0] # first character
        return (codon_nt == "A" and trna_nt in ['A','I','M','?']) or \
                (codon_nt == "C" and trna_nt in ['A','U','P','I','?','Q']) or \
                (codon_nt == "G" and trna_nt in ['A','U','&','3','1','~','N','S',')','{','V','P','?','M']) or \
                (codon_nt == "U" and trna_nt in ['G','#','W','U','V','P','I','Q'])

    cognate_wc_matrix = np.zeros((len(t_rnas.anticodon), len(codons.codon)))
    cognate_wobble_matrix = np.zeros((len(t_rnas.anticodon), len(codons.codon)))
    nearcognate_matrix = np.zeros((len(t_rnas.anticodon), len(codons.codon)))

    if verbose:
        print("Populating WC matrix...")
    # populate cognate WC matrix if WC criteria matched
    for anticodon_index, anticodon in enumerate(t_rnas.anticodon):
        for codon_index, codon in enumerate(codons.codon):
            if second_wc(codon, anticodon) and first_wc(codon, anticodon) and third_wc(codon, anticodon):
                cognate_wc_matrix[anticodon_index, codon_index] = 1
    if verbose:
        print("done.")
        print("Populating wobble matrix...")


    #populate cognate wobble matrix if wobble criteria matched, amino acid is correct, and WC matrix entry is 0
    #if incorrect amino acid, assign to near-cognates

    for anticodon_index, anticodon in enumerate(t_rnas.anticodon):
        for codon_index, codon in enumerate(codons.codon):
            if cognate_wc_matrix[anticodon_index,codon_index] == 0 and second_wc(codon, anticodon) and\
               first_wc(codon, anticodon) and third_wobble(codon, anticodon):
                if t_rnas["three.letter"][anticodon_index] == codons["three.letter"][codon_index]:
                    cognate_wobble_matrix[anticodon_index,codon_index] = 1
                else:
                    nearcognate_matrix[anticodon_index,codon_index] = 1

    if verbose:
        print('done.')
        print('Populating nearcognate matrix...')

    #populate near-cognate matrix if:
    #wobble and WC matrix entries are 0,
    #wobble criteria are matched

    for anticodon_index, _ in enumerate(t_rnas.anticodon):
        for codon_index, _ in enumerate(codons.codon):
            if (cognate_wc_matrix[anticodon_index,codon_index] == 0 and\
                cognate_wobble_matrix[anticodon_index,codon_index] == 0) and\
                (second_wc(codons.codon[codon_index],t_rnas.anticodon[anticodon_index]) and \
                 (first_wobble(codons.codon[codon_index],t_rnas.anticodon[anticodon_index]) or\
                     first_wc(codons.codon[codon_index],t_rnas.anticodon[anticodon_index])) and\
                 (third_wobble(codons.codon[codon_index],t_rnas.anticodon[anticodon_index]) or\
                     third_wc(codons.codon[codon_index],t_rnas.anticodon[anticodon_index])) ):
                nearcognate_matrix[anticodon_index,codon_index] = 1

    if verbose:
        print('done.')

    #Sanity checks

    #Check whether any tRNA:codon combination is assigned 1 in more than one table (this should not occur)

    testsum = cognate_wc_matrix + cognate_wobble_matrix + nearcognate_matrix
    if np.any(testsum>1):
        print('Warning: multiple relationships for identical tRNA:codon pairs detected.')
        return {}
    elif verbose:
        print('No threesome errors detected.')

    return {"cognate.wc.matrix":cognate_wc_matrix, "cognate.wobble.matrix":cognate_wobble_matrix,
            "nearcognate.matrix":nearcognate_matrix}

def plot_matrix(matrices_dict, t_rnas, codons, save_fig = None):
    """
    Plots the pairing matrices.
    """

    colours=['g', 'y', 'r']
    labels = list(matrices_dict.keys())
    i = 0
    plt.figure(figsize=(25,15))
    plt.grid(True)
    for k in matrices_dict.keys():
        matches_dict = np.argwhere(matrices_dict[k] == 1)#[:,1]
        # display(c)
        plt.plot(matches_dict[:,1], matches_dict[:,0], colours[i] + 's', label=labels[i])
        i +=1
    plt.xticks(range(len(codons.codon)), codons.codon, rotation = 45)
    plt.yticks(range(len(t_rnas.anticodon)), t_rnas.anticodon)
    plt.legend()
    if save_fig is not None:
        plt.savefig(save_fig)
    plt.show()


def make_concentrations(matrices_dict, t_rnas, codons, concentration_col_name = 'gene.copy.number', total_trna=190):
    """
    Given a tRNA matrix, and the decoding matrix, generates a concentrations DataFrame.

    TRNAs: DataFrame with the tRNAs: anticodon, gene.copy.number
    Matrices: pairing matrices generated by make_matrix
    Codons: DataFrame with the codons to be used
    concentration_col_name: name of the concentrations column name. Default = 'gene.copy.number'
    total_Trna: default value is 190 (micro moles).
    """
    wc_cognate = matrices_dict["cognate.wc.matrix"]
    wobblecognate = matrices_dict["cognate.wobble.matrix"]
    nearcognate = matrices_dict["nearcognate.matrix"]

    # construct empty results dataframe
    trna_concentrations = pd.DataFrame(codons[[codons.columns[0],codons.columns[2]]])
    trna_concentrations["WCcognate.conc"] = 0.0
    trna_concentrations["wobblecognate.conc"] = 0.0
    trna_concentrations["nearcognate.conc"] = 0.0

    #calculate a conversion factor to convert the abundance factor to a molar concentration
    print('using: '+ concentration_col_name)
    conversion_factor = np.float64(total_trna / np.float64(t_rnas[concentration_col_name].sum()) * 1e-6)

    # go through the WCcognates matrix and for each entry of 1 add the abundance of the tRNA from the abundance table
    # to the concentration table
    for codon_index, _ in enumerate(codons.codon):
        for anticodon_index, _ in enumerate(t_rnas.anticodon):
            if wc_cognate[anticodon_index, codon_index] == 1:
                trna_concentrations.loc[codon_index, "WCcognate.conc"] =\
                    trna_concentrations["WCcognate.conc"][codon_index] +\
                        (t_rnas[concentration_col_name][anticodon_index]*conversion_factor)

    #ditto for wobblecognate
    for codon_index, _ in enumerate(codons.codon):
        for anticodon_index, _ in enumerate(t_rnas.anticodon):
            if wobblecognate[anticodon_index, codon_index] == 1:
                trna_concentrations.loc[codon_index, "wobblecognate.conc"] =\
                    trna_concentrations["wobblecognate.conc"][codon_index] +\
                        (t_rnas[concentration_col_name][anticodon_index]*conversion_factor)

    #ditto for nearcognates
    for codon_index, _ in enumerate(codons.codon):
        for anticodon_index, _ in enumerate(t_rnas.anticodon):
            if nearcognate[anticodon_index, codon_index] == 1:
                trna_concentrations.loc[codon_index, "nearcognate.conc"] =\
                    trna_concentrations["nearcognate.conc"][codon_index] +\
                        (t_rnas[concentration_col_name][anticodon_index]*conversion_factor)
    return trna_concentrations

##example of how to use:

# tRNAs = pd.read_csv('/home/heday/Projects/R3/Native Spike and B117 Kent/HEK293_processed.csv')
# codons = pd.read_csv('/home/heday/Projects/R3/Native Spike and B117 Kent/codons.csv')
# matrix = make_matrix(tRNAs, codons, verbose=True)
