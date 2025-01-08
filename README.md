# Elongation Simulator

This is the repository containing the source code of the Protein Synthesis.

It contains python scripts to calculate tRNA concentrations and two simulators:
+ codon_simulator
+ sequence_simulator

## Codon Simulator
This class is essentially an implementation of the Gillespie algorithm, and allows running stochastic simulations of the decoding of individual codons efficiently.

## Sequence Simulator
This class relies on a modified implementation of the Gillespie algorithm. This simulator tracks ribosome positional information, allowing the simulation of mRNA transcripts containing any number of elongating ribosomes. It also allows for setting their initiation and termination rates, and a choice of criteria to stop the simulations.

