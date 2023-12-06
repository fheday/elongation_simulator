#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 02 11:32:41 2019

Re-implementation of the R-code that calculates concentrations
This script is intented to generate new concentrations (WC, wobble, near) with different
interpretation of wobble and near cognates. it is intended to be use in GA algorithm.

@author: heday
"""

from importlib import resources
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import concentrations

def identify_3_base_pairing(codon: str, trna: str, basepairing_rules:dict) -> [str]:
    result = []
    for i, c in enumerate(codon):
        anticodon = trna[2 - i]
        if anticodon in basepairing_rules["Watson-Crick"][c]:
            result.append("WC")
        elif anticodon in basepairing_rules["Wobble"][c]:
            result.append("Wo")
        else:
            result.append("X")
    return result

def is_near_cognate(pairing_id:[str], near_cognate_rules:[[str]]) -> bool:
    for partial_rule in near_cognate_rules:
        matches = 0
        for i, pairing_type in enumerate(pairing_id):
            if partial_rule[i].upper() == 'X':
                # doesn't matter, it is always a match.
                matches += 1
            elif pairing_type.upper() == partial_rule[i].upper():
                # match.
                matches += 1
            elif pairing_type.upper() == 'WC' and partial_rule[i].upper() == 'WO':
                # this is also a match.
                matches += 1
        if matches == 3:
            return True
    return False


def make_matrix(t_rnas: pd.DataFrame, codons: pd.DataFrame, verbose=False, settings_file_name="default_basepairing.json"):
    """
    Given a DataFrame with tRNA concentrations, and another DataFrame with codons information,
    generates a decoding matrix

    TRNAs: DataFrame with the tRNAs: anticodon, gene.copy.number
    Codons: DataFrame with the codons to be used
    verbose: print information about the building stages as it goes
    settings_file: json file with the basepairing rules.
    """
    # open file with basepairing rules
    if settings_file_name == "default_basepairing.json":
        with resources.path(concentrations, settings_file_name) as re:
            with open(re, "r") as file:
                basepairing_rules = json.load(file)
    else:
        with open(settings_file_name, "r") as file:
                basepairing_rules = json.load(file)
    # sanity checks
    if "Watson-Crick" not in basepairing_rules.keys():
        raise ValueError("settings file does not contain Watson-Crick pariring rules")
    elif "Wobble" not in basepairing_rules.keys():
        raise ValueError("settings file does not contain Wobble pariring rules")
    elif "A" not in basepairing_rules["Watson-Crick"]:
        raise ValueError("settings file does not rules for Watson-Crick A basepairing")
    elif "C" not in basepairing_rules["Watson-Crick"]:
        raise ValueError("settings file does not rules for Watson-Crick C basepairing")
    elif "G" not in basepairing_rules["Watson-Crick"]:
        raise ValueError("settings file does not rules for Watson-Crick G basepairing")
    elif "U" not in basepairing_rules["Watson-Crick"]:
        raise ValueError("settings file does not rules for Watson-Crick U basepairing")
    elif "A" not in basepairing_rules["Wobble"]:
        raise ValueError("settings file does not rules for Wobble A basepairing")
    elif "C" not in basepairing_rules["Wobble"]:
        raise ValueError("settings file does not rules for Wobble C basepairing")
    elif "G" not in basepairing_rules["Wobble"]:
        raise ValueError("settings file does not rules for Wobble G basepairing")
    elif "U" not in basepairing_rules["Wobble"]:
        raise ValueError("settings file does not rules for Wobble U basepairing")
    # settings seems fine. Proceed.

    # check if tRNAs have 'anticodon' column
    if 'anticodon' not in t_rnas.columns:
        print('tRNA list must contain a column named "anticodon".')
        return
    # check if codons have 'codon' column
    if 'codon' not in codons.columns:
        print('Codon list must contain a column named "codon".')
        return



    cognate_wc_matrix = np.zeros((len(t_rnas.anticodon), len(codons.codon)))
    cognate_wobble_matrix = np.zeros((len(t_rnas.anticodon), len(codons.codon)))
    nearcognate_matrix = np.zeros((len(t_rnas.anticodon), len(codons.codon)))

    if verbose:
        print("Populating WC matrix...")
    # populate cognate WC matrix if WC criteria matched
    for anticodon_index, anticodon in enumerate(t_rnas.anticodon):
        for codon_index, codon in enumerate(codons.codon):
            if identify_3_base_pairing(codon, anticodon, basepairing_rules) == ["WC", "WC", "WC"]:
                cognate_wc_matrix[anticodon_index, codon_index] = 1

    if verbose:
        print("done.")
        print("Populating wobble matrix...")


    #populate cognate wobble matrix if wobble criteria matched, amino acid is correct, and WC matrix entry is 0
    #if incorrect amino acid, assign to near-cognates
    for anticodon_index, anticodon in enumerate(t_rnas.anticodon):
        for codon_index, codon in enumerate(codons.codon):
            if cognate_wc_matrix[anticodon_index,codon_index] == 0 and\
               identify_3_base_pairing(codon, anticodon, basepairing_rules) == ["WC", "WC", "Wo"]:
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
                is_near_cognate(identify_3_base_pairing(codons.codon[codon_index], t_rnas.anticodon[anticodon_index], basepairing_rules), basepairing_rules["Pairing Rules"]["Near-Cognate"]["base-level"]):
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

def plot_matrix(matrices_dict: dict, t_rnas: pd.DataFrame, codons: pd.DataFrame, save_fig = None):
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

def print_codon_anticodon_pairings(matrices_dict: dict, t_rnas: pd.DataFrame, codons: pd.DataFrame):
    """
    prints the pairing of all codons with all anticodons.
    """
    for k in matrices_dict.keys():
        print(" Pairings for: " + k)
        matches_dict = np.argwhere(matrices_dict[k] == 1)
        for anticodon, c in matches_dict:
            print (codons.codon[c] + " - " + t_rnas.anticodon[anticodon])

def make_concentrations(matrices_dict: dict, t_rnas: pd.DataFrame, codons: pd.DataFrame, concentration_col_name = 'gene.copy.number', total_trna=190, verbose=False):
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
    trna_concentrations = pd.DataFrame(codons[[codons.columns[0]]])
    trna_concentrations["WCcognate.conc"] = 0.0
    trna_concentrations["wobblecognate.conc"] = 0.0
    trna_concentrations["nearcognate.conc"] = 0.0

    #calculate a conversion factor to convert the abundance factor to a molar concentration
    if verbose:
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
