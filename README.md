
Explore dN/dS ratios in viral populations simulated with SANTA

The idea is to verify that selection functions in SANTA are working by
measuring dN/dS ratios.

dNdS ratios are calculated by PAML for sequences arrange on a tree. It is
for that reason that we use BEAS?T to infer a phylogenetic tree on the
simulated sequences.

`scons` to generate all the data.
$ scons  build/2.5E-5/frequency/indel/pop:1000/rep:9/gen:5000/N:10/strict

First we ru santa for a "long" time to generate a diverse population
of sequences sampled at various timepoints.
Then choose two timepoints at which we will extract samples to work
with.
Use `mafft` to align the samples because some of our simulations
include indels.
`beast` infers a phylogenetic tree for the combined sequences.
`treeannotator` calculates an MCC tree.

After that we can start running PAML to calculate dN/dS ratios.

$ bin/paml.py -vv codeml build/2.5E-5/frequency/noindel/pop:1000/rep:9/gen:5000/N:10/sample.fa build/2.5E-5/frequency/noindel/pop:1000/rep:9/gen:5000/N:10/strict/mcc.tree


