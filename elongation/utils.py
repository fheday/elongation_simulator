'''
Utility functions
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
from collections import Counter

"""
For each entry in the simulation log, returns the number of elongating ribosomes.
"""
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
