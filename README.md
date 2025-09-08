# Elongation Simulator

Elongation Simulator is package for simulating Protein Synthesis.

It contains python scripts to calculate tRNA concentrations and two simulators:
+ codon_simulator
+ sequence_simulator

## Installation from PyPi Binary Package


```python
pip install elongation-simulator
```

## Installation from Source

The source installation works on Linux, MacOS (Intel and Apple). It needs a C++17 compatible compiler. The compilation is done through pip:
```python
pip install .
```

## Codon Simulator
This class is essentially an implementation of the Gillespie algorithm, and allows running stochastic simulations of the decoding of individual codons efficiently. A simple use can be find below:

```python
from codon_simulator import CodonSimulator
sim = CodonSimulator()
esim.load_concentrations(sim.saccharomyces_cerevisiae_concentrations) # the simulator already comes with the yeast concentrations
sim.set_codon_for_simulation('AAG') # sets the codon to be simulated
sim.run_repeatedly_get_average_time(1000) # simulates the codon 1k times and returns its average decoding time
```

## Sequence Simulator
This class relies on a modified implementation of the Gillespie algorithm. This simulator tracks ribosome positional information, allowing the simulation of mRNA transcripts containing any number of elongating ribosomes. It also allows for setting their initiation and termination rates, and a choice of criteria to stop the simulations. Below is an example of a simple use of the Sequence Simulator:

```python
from sequence_simulator import SequenceSimulator
sim = SequenceSimulator()
sim.load_concentrations(sim.saccharomyces_cerevisiae_concentrations) # the simulator already comes with the yeast concentrations
sim.input_MRNA("ATGTTCAGCGAATTAATTAACTTCCAAAATGAAGGTCATGAGTGCCAATGCCAATGTGGTAGCTGCAAAAATAATGAACAATGCCAAAAATCATGTAGCTGCCCAACGGGGTGTAACAGCGACGACAAATGCCCCTGCGGTAACAAGTCTGAAGAAACCTGA")
sim.set_initiation_rate(100) # sets the rate of initiation to 100 ribosomes per second.
sim.set_termination_rate(100) # sets the rate for ribosomes terminations to 100 terminations per second.
sim.set_finished_ribosomes(2) # stops the simulation after 2 ribosomes successfully terminated.
sim.run() # runs the simulation
sim.dt_history # displays the times where ribosomes moved: either initiation, elongation or termination.
sim.ribosome_positions_history # displays the postiions of all elongating ribosomes during the simulation.
```
