# SANTA-dNdS

Explore dN/dS ratios in viral populations simulated with SANTA.
The idea is to verify that selection functions in SANTA are working
through standard measures for purifying and
positive selection.

## Basic usage

To gather statistics and build plots that explore [SANTA](https://github.com/santa-dev/santa-sim) memory and runtime performance,
```
	$ scons -j 10
```
Warning: No attempt has been made to make this code portable to
environments outside the Hutch.

Running the simulations under various combinations of population size
and sequence length will take many hours to complete.  Two
optimizations are included in the `SConstruct` file to hopefully make this process faster.  The
`SConstruct` file uses the `srun` command from the
`[SLURM](http://slurm.schedmd.com/)` cluster management package to launch jobs on the Hutch center cluster resources.  Also,
because the heavy lifting is offloaded to the cluster, we use `-j N`
to launch multiple simulations in parallel.

Configuration files and output from the inividual simulations are left
directories under `build/`.  Aggregated statistics are left in
`build/results.csv`.  Those results can be plotted with
`plot_results.R`

# Results

![Fitness plots](figures/fitness_plots.png)

# Methodology

1. Sample from simulated populations under a variety of conditions.
2. Measure dN/dS ratios to  determine if simulation conditions affect
sample populations in predictable ways.

We configure SANTA to evolve populations under four selection
conditions: 'noselection', 'purifyingchem', 'empiricalvalues',
'empiricalvalues_homoresidue'.  All of these conditions are unbiased in
the generation of substitution mutations, and all configurations
penalize in-frame STOP codons as lethal mutations.  Further details of
each configuration are explained in the next section.

Multiple replicated simulations
are configured and run for each of three fitness values: low, medium,
and high.  The actual fitness values
depend upon the particular selection condition being use, i.e. low
might be '0.001' for 'purifyingchem' and '.1' for  'empiricalvalues'.
The simulations are configured to emit samples of sequences at multiple
intervals to assess how the population changes over time.

The sequences sampled from each simulation are deduplicated before
being used as input to [BEAST](http://beast.bio.ed.ac.uk/) to infer a
phylogenetic tree.  The tree and deduplicated samples are then passed
to [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) for
calculation of a single dN/dS value for the entire tree.

# Details of Selection configurations

'noselection' has no fitness constraints other than
penalizing in-frame STOP codons.

'purifyingchem' uses an `<purifyingFitness>` element to apply a linear
decay fitness function at each residue position which favors residues
that are chemically similar to the starting residue.
```xml
	 <fitnessFunction>
	     <purifyingFitness>
	         <feature>C2C3</feature>
	         <sites>1-156</sites>
	         <rank>
	             <order>chemical</order>
	             <breakTies>random</breakTies>
	         </rank>
	         <fitness>
	             <lowFitness>$fitness</lowFitness>
	             <minimumFitness>0.1</minimumFitness>
	         </fitness>
	     </purifyingFitness>
	 </fitnessFunction>
```

In other words,
each residue will have a preference for the original residue at that
site, with a decaying preference for chemically similar residues.
Outside of chemically similar residues, fitness drops off to a constant,
minimal value.  This should serve to maintain the amino acids content
of a simulated genome, with some variation due to  synonymous
mutations and chemically similar
residues.

'empiricalvalues' uses an `<empiricalSelection>` element that assigns
specific values to each amino acid residue.  For this experiment all
residues except for Lysine get asssigned the same fitness value.  the
fitness value assigned to Lysine is a parameter of the simulation.

```xml
	 <fitnessFunction>
	   <empiricalFitness>
	     <feature>C2C3</feature>
	     <sites>1-156</sites>
	     <!-- assign fitness to each amino acid at this site -->
	     <!-- Order of Amino cacids is A-C-D-E-F-G-H-I-K-L-M-N-P-Q-R-S-T-V-W-Y -->
	     <values>
	 	  0.010 0.010 0.010 0.010 0.010
	 	  0.010 0.010 0.010 $fitness 0.010
	 	  0.010 0.010 0.010 0.010 0.010
	 	  0.010 0.010 0.010 0.010 0.010
	 	</values>
	   </empiricalFitness>
	 </fitnessFunction>
```

'empiricalvalues_homoresidue' uses the same fitness
function as'empiricalvalues', but the inital population starts with a
genomethat entrely consisting of codons that code for Lysine, e.g.
```xml
    <sequences>
	AAAAAAAAAAAAAAAAAAAAAAAA...AAAAAAAAAAAAAAAAAAAAAAA
    </sequences>

```



