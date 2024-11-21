'''
Utility functions
'''

import numpy as np

def get_number_elongating_ribosomes_per_time(sim):
    """
    For each entry in the simulation log, returns the number of elongating ribosomes.
    """
    return [len(entry) for entry in sim.ribosome_positions_history], np.cumsum(sim.dt_history)


def get_codon_average_occupancy(sim):
    '''
    Returns a tuple where the first element is a vector with the
    enlogation duration of the ribosomes that terminated in the simulation, and
    the second element is a vector with the iteration where such ribosomes
    started enlogating.
    '''
    # sim.updateRibosomeHistory()
    return sim.getEnlogationDuration()
