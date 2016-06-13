# diblockcopoly
=======================================

This is a function that uses polymer field theory to find phase behavior of diblock copolymers

The polymers are modeled as wormlike chains, Gaussian chains, and perfectly rigid rods.
Phase transition spinodal and critical wavemode of phase segregation can be found at different
chemical correlation and monomer rigidities.

Renormalized phase diagrams are found by F-H/Brazovskii theory of free energy expansion
up to quartic-order density fluctuations.

## Run example
example.m provides simple demo that utilizes the functions including
- plotphase, plots a mean-field diblock copolymer phase diagram
- plotphaseRG, plots a diblock copolymer phase diagram with density fluctuations
- spinodal, calculates mean-field phase transition spinodal of diblock copolymer
- spinodalRG, calculates renormalized phase transition spinodal of diblock copolymer with density flucutations
- densityRG, plots density-density correlations with mean-field and flucutation theory
- calcgamma, calculates vertex functions of free energy expansion

## Input specifications
N, number of Kuhn steps of total chain. In general chains are modeled as worm-like chains. In the limit N>1e4, Gaussian chain statistics are used; in the other limit N<1e-4, perfectly rigid rod statistics are used
FA, chemical composition of A type monomers