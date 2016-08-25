#!/usr/bin/env scons
# -*- coding: utf-8 -*-

'''
Measure the performance of founder inference various numbers of sequences.

This SCons file varies clock model (strict or relaxed)

I think what we learned form this simulation is that static vs relaxed
clock doesn't make much difference.  Literature and other sources
suggest using strict clock for intraspecific samples, which is what we
have here.  We will move forward with just the strict clock.


'''



import os
import sys
import os.path
import numpy as np

from nestly import Nest
from nestly.scons import SConsWrap
from SCons.Script import Environment

from sconsutils import Wait
import sconsutils
import subprocess

environ = os.environ.copy()

env = Environment(ENV=environ)
env.PrependENVPath('PATH', 'bin')
env['PRANK']='/home/cwarth/src/matsen/prank/src/prank'
env['SANTAJAR']= os.path.expanduser('~matsengrp/local/lib/santa.jar')
env['SANTAJAR']= os.path.expanduser('~/src/matsen/santa-wercker/dist/santa.jar')

env['LONGEVITY'] = 20000	# number of generations to run SANTA simulation.

def cmd_exists(cmd):
    return subprocess.call("type " + cmd, shell=True, 
        stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

# check if `BEAST` is available, but only if we aren't in a dry-run.
if not GetOption('no_exec') and not cmd_exists('beast'):
    msg = '''
    `beast` not found on PATH
     Consider using,
        module load BEAST/1.8.2
    '''
    print(msg)
    sys.exit(1)

n = Nest(base_dict={})
w = SConsWrap(n, 'build', alias_environment=env)

# adding aggregate for running build_graph.py
w.add_aggregate('resultsList', list)

w.add('mutationrate', ['2.5E-5'], create_dir=False)

w.add('indel_model', ['noindel'], create_dir=False)

w.add('nseqs', [10], label_func=lambda n: 'N_'+str(n), create_dir=False)

w.add('population', [1000], label_func=lambda p: 'pop_'+str(p), create_dir=False)

w.add('selection_model', ['empiricalvalues_homoresidue', 'noselection', 'purifyingchem', 'empiricalvalues']) 

def fitness_range(c): # takes control dictionary
    # choose selection parameters based on selection_model value.
    # some selection models require a different range of fitness values.
    # some also vary the starting sequence.

    return {
        'noselection': np.linspace(0.0100, 0.02,3),
        'purifyingchem': np.linspace(0.0100, 1, 3),
        'empiricalvalues': np.linspace(0.0100, 0.02, 3),
        'empiricalvalues_homoresidue': np.linspace(0.0100, 0.0105, 3)
    }[c['selection_model']]


w.add('fitness', fitness_range, label_func=lambda n: 'fit_'+str(n))

def sequence_files(c): # takes control dictionary
    return {
        'noselection': ['templates/HIV1C2C3.fasta'],
        'purifyingchem': ['templates/HIV1C2C3.fasta'],
        'empiricalvalues': ['templates/HIV1C2C3.fasta'],
        'empiricalvalues_homoresidue': ['templates/homoresidue.fasta']
    }[c['selection_model']]

w.add('sequence', sequence_files, create_dir=False)

def template_file(c): # takes control dictionary
    return {
        'noselection': ['templates/santa_${selection_model}_${indel_model}.template'],
        'purifyingchem': ['templates/santa_${selection_model}_${indel_model}.template'],
        'empiricalvalues': ['templates/santa_${selection_model}_${indel_model}.template'],
        'empiricalvalues_homoresidue': ['templates/santa_empiricalvalues_noindel.template']
    }[c['selection_model']]

w.add('template', template_file, create_dir=False)


@w.add_target_with_env(env)
def santa_config(env, outdir, c):
    return env.Command(
        os.path.join(outdir, "santa_config.xml"),
		['${template}', '${sequence}'],
		"mksanta.py  -p patient1 ${SOURCES}   >${TARGET}"
        )[0]



w.add('replicates', range(5), label_func=lambda r: 'rep_'+str(r))

@w.add_target_with_env(env)
def santa_lineage(env, outdir, c):
    return env.Command(os.path.join(outdir, "donorlineage.fa"),
                       [ c['santa_config'], env['SANTAJAR'] ],
                       [  # santa will produce output files in its current directory.
                          # so need to change to output directory before execution.
                          Copy('${OUTDIR}/santa_config.xml', '${SOURCES[0]}'),
                          'cd ${OUTDIR} && srun java -mx512m -jar ${SOURCES[1]} -mutationrate=${mutationrate} -population=${population} -longevity=${LONGEVITY} -fitness=${fitness} santa_config.xml',
                          Copy('${TARGET}', '${OUTDIR}/santa_out.fa')
                       ])[0]

w.add('timepoint', [500, 5000, 20000], label_func=lambda p: 'gen_'+str(p))

## Extract the founder sequence from the santa config file into a FASTA file.
## This makes it easier for the distance.py script to grab it for comparison.
@w.add_target_with_env(env)
def sample(env, outdir, c):
    return env.Command(
        os.path.join(outdir, 'sample.fa'),
        [ c['santa_lineage'] ],
        [
            # Append fake dates to the sequence ids..
            # The sequence ids will be parsed when building the beast config file and tip dates will be created to matched the dates on the sequences.
			'''
			fasta_sample.py --fasta-file ${SOURCES[0]} --n-sequences ${nseqs} --pattern '_${timepoint}_' | \
			seqmagick convert --pattern-replace '^([^\|]*)\|.*$' '\\1|1M|XXX|XXX|2011_11_10' - - >${TARGET}
			'''
        ])[0]


# # align sample
# @w.add_target_with_env(env)
# def align(env, outdir, c):
#     return env.Command(
#         os.path.join(outdir, 'sample_aln.fa'),
#         [ c['sample'] ],
#         'mafft --quiet --auto ${SOURCE} >${TARGET}')[0]

@w.add_target_with_env(env)
def fasta2phylip(env, outdir, c):
    return env.Command( os.path.join(outdir, 'sample_aln.phylip'),
                        c['sample'],
                        'fasta2phylip.py ${SOURCE} ${TARGET}'
                        )[0]

# create the BEAST config file from sequences extracted from two patient simulations
@w.add_target_with_env(env)
def config_beast(env, outdir, c):
    return env.Command(os.path.join(outdir, 'beast_in.xml'),
                       [ 'templates/beast_strict.template', c['sample']],
                       "mkbeast_rv217.py  --template  ${SOURCES[0]} ${SOURCES[1]}  >${TARGET}")[0]


@w.add_target_with_env(env)
def runbeast(env, outdir, c):
    target = [ os.path.join(outdir, 'beastout.log'),
               os.path.join(outdir, 'beastout.trees'),
               os.path.join(outdir, 'beastcmd.log'),
               os.path.join(outdir, 'srun.log')
                   ]
    return env.Command(target,
                       c['config_beast'],
                       [ "srun --time=30 --chdir=${OUTDIR} --output=${TARGETS[2]} beast -overwrite -beagle ${SOURCE.file} >${TARGETS[3]} 2>&1",
                         Wait(target)
                       ])[1]

@w.add_target_with_env(env)
def mcc(env, outdir, c):
    return env.Command(os.path.join(outdir, 'mcc.nexus'),
                        c['runbeast'],
                       'treeannotator ${SOURCE} >${TARGET} ')

@w.add_target_with_env(env)
def nexus2newick(env, outdir, c):
    return env.Command(os.path.join(outdir, 'mcc.newick'),
                        c['mcc'],
                       'nexus2newick.py ${SOURCES} ${TARGET} ')

@w.add_target_with_env(env)
def config(env, outdir, c):
    return env.Command(os.path.join(outdir, 'codeml.ctl'),
                        [  'codeml.ctl', c['nexus2newick'], c['fasta2phylip']],
                        [
                            "sed 's#seqfile.txt#${SOURCES[2].file}#g' <${SOURCES[0]} >${TARGET}"
                        ])[0]

    
@w.add_target_with_env(env)
def results(env, outdir, c):
    return env.Command(os.path.join(outdir, 'results.txt'),
                       c['config'],
                       'cd ${SOURCE.dir} && codeml ${SOURCE.file}'
                       )[0]

@w.add_target_with_env(env)
def aggregate(env, outdir, c):
    c['resultsList'].append(c['results'])

w.pop('mutationrate')

@w.add_target_with_env(env)
def collect(env, outdir, c):
    return env.Command(os.path.join(outdir, 'results.csv'),
                       c['resultsList'],
                       'parseresults.py -o ${TARGET} ${SOURCES}'
                       )




