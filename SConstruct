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
import os.path
import glob
import tempfile
import time

from nestly import Nest
from nestly.scons import SConsWrap
from SCons.Script import Environment
from SCons.Action import ActionFactory
import SCons.Util


# Running commands on the cluster sometimes has the unfortunate side-effect of
# letting distributed filesystems get out of sync.  A file that is written on
# the cluster may not be visible on local machines for several seconds.  This
# doesn't happen all the time, and it can be difficult to demonstrate, but it
# often occurs when local and cluster commands manipulate the same file in
# quick succession.  The Wait() action is meant to wait for a file to become
# locally visible after running a command on the cluster (via srun or salloc).
#
# Wait() usage is typically,
# 	target='output.file'
# 	env.Command(target, 'source.file',
#   	        [ "srun some-command <${TARGET}",
#    			   Wait(target)
#				])
#
# This will cause the execution to pause after running 'some-command' until the target shows up on the local machine.
# The target will be polled on a 2-second interval, and the command will fail if the target does not show up within about 10 seconds.




def get_paths_str(dest):
    # If dest is a list, we need to manually call str() on each element
    if SCons.Util.is_List(dest):
        elem_strs = []
        for element in dest:
            elem_strs.append('"' + str(element) + '"')
        return '[' + ', '.join(elem_strs) + ']'
    else:
        return '"' + str(dest) + '"'

# https://github.com/azatoth/scons/blob/73f996e59902d03ec432cc662252aac5fb72f1f8/src/engine/SCons/Defaults.py 
def wait_func(dest):
    SCons.Node.FS.invalidate_node_memos(dest)
    if not SCons.Util.is_List(dest):
        dest = [dest]
    for entry in dest:
        count = 0
        limit = 20
        while not os.path.isfile(entry) or os.stat(entry).st_size == 0:
            print("waiting for {}...".format(entry))
            time.sleep(2)
            count = count + 1
            if count > limit:
                print("failing wait for {}...".format(entry))
                return 1
    return 0

Wait = ActionFactory(wait_func, lambda dir: 'Wait(%s)' % get_paths_str(dir))

environ = os.environ.copy()

env = Environment(ENV=environ)
env.PrependENVPath('PATH', '../bin')
env['PRANK']='/home/cwarth/src/matsen/prank/src/prank'
env['SANTAJAR']= os.path.expanduser('~matsengrp/local/lib/santa.jar')
env['SANTAJAR']= os.path.expanduser('~/src/matsen/santa-dev/dist/santa.jar')

env['LONGEVITY'] = 25000	# number of generations to run SANTA simulation.

n = Nest(base_dict={})
w = SConsWrap(n, 'build', alias_environment=env)


n.add('mutationrate', ['2.5E-5'])

n.add('selection_model', ['noselection', 'purifying',  'frequency'])
n.add('indel_model', ['noindel', 'indel'])


@w.add_target_with_env(env)
def santa_config(env, outdir, c):
    return env.Command(os.path.join(outdir, "santa_config.xml"),
                       ['templates/santa_{selection_model}_{indel_model}.template'.format(**c), 'templates/HIV1C2C3.fasta'],
                       "mksanta.py  -p patient1 ${SOURCES}   >${TARGET}")[0]


n.add('population', [1000], label_func=lambda p: 'pop='+str(p))

n.add('replicates', range(10), label_func=lambda r: 'rep='+str(r))

@w.add_target_with_env(env)
def santa_lineage(env, outdir, c):
    return env.Command(os.path.join(outdir, "donorlineage.fa"),
                       [ c['santa_config'], env['SANTAJAR'] ],
                       [  # santa will produce output files in its current directory.
                          # so need to change to output directory before execution.
                          Copy('${OUTDIR}/santa_config.xml', '${SOURCES[0]}'),
                          'cd ${OUTDIR} && srun --time=30 --output=srun.log java -mx512m -jar ${SOURCES[1]} -mutationrate=${mutationrate} -population=${population} -longevity=${LONGEVITY} santa_config.xml',
                          Copy('${TARGET}', '${OUTDIR}/santa_out.fa')
                       ])[0]

n.add('timepoint', [1000, 5000], label_func=lambda p: 'gen='+str(p))

n.add('nseqs', [10], label_func=lambda n: 'N='+str(n))

## Extract the founder sequence from the santa config file into a FASTA file.
## This makes it easier for the distance.py script to grab it for comparison.
@w.add_target_with_env(env)
def sample(env, outdir, c):
    target = os.path.join(outdir, 'sample.fa'.format(**c))
    fasta_sample = r'fasta_sample.py --fasta-file ${{SOURCES[0]}} --n-sequences {} --pattern "_{}_"'
    mogrify = r'seqmagick convert --pattern-replace "^([^\|]*)\|.*$" "\1|{}" - -'
    samplecmd = fasta_sample + '|' + mogrify 
    return env.Command(target,
                [ c['santa_lineage'] ],
                [
                    samplecmd.format(c['nseqs'], c['timepoint'], '1M|XXX|XXX|2011/11/10') + ' >${TARGET}',
                    samplecmd.format(c['nseqs'], str(c['timepoint']+500), '6M|XXX|XXX|2012/03/02') + ' >>${TARGET}',
#                    r'fasta_sample.py --fasta-file ${SOURCES[0]} --n-sequences 10 --pattern "_${timepoint}_" | seqmagick convert --pattern-replace "^([^\|]*)\|.*$" "\1|1M|05WG|NFLG|2011/11/10" - - >${TARGET}',
#                    r'fasta_sample.py --fasta-file ${SOURCES[0]} --n-sequences 12 --pattern "_' + str(c['timepoint']+500) + r'_" | seqmagick convert --pattern-replace "^([^\|]*)\|.*$" "\1|6M|08RH|RH|2012/03/02" - - >>${TARGET}'
                ])[0]


# align sample
@w.add_target_with_env(env)
def multiple_alignment(env, outdir, c):
    #founder = os.path.join(outdir, 'sample_dedup_aln.fa'.format(**c))
    target = '{}_aln.fa'.format(os.path.splitext(str(c['sample']))[0])
    cmd = 'mafft --quiet --auto ${SOURCE} >${TARGET}'

    return env.Command(target,
                [ c['sample'] ],
                [ cmd ])


# For each data set, we evaluated both strict and relaxed uncorrelated lognormal clock models.
# These are captured in the different template files used to build the BEAST configs.
#
# See: McCloskey, R. M., Liang, R. H., Harrigan, P. R., Brumme, Z. L.,
# & Poon, A. F. Y. (2014). An Evaluation of Phylogenetic Methods for
# Reconstructing Transmitted HIV Variants using Longitudinal Clonal
# HIV Sequence Data. Journal of Virology, 88(11),
# 6181–94. doi:10.1128/JVI.00483-14

n.add('clock_model', ['relaxed', 'strict'])	


# create the BEAST config file from sequences extracted from two patient simulations
@w.add_target_with_env(env)
def config_beast(env, outdir, c):
    target = os.path.join(outdir, 'beast_in.xml')
        
    cmd = ("mkbeast_rv217.py  --template  ${SOURCES[0]} ${SOURCES[1]}  >${TARGET}")

    return env.Command(target,
                       [ '../templates/beast_{}.template'.format(c['clock_model']), c['multiple_alignment']],
                       cmd)[0]


@w.add_target_with_env(env)
def runbeast(env, outdir, c):
    target = [ os.path.join(outdir, 'ancestralSequences.log'),
               os.path.join(outdir, 'beastout.log'),
               os.path.join(outdir, 'beastout.trees'),
               os.path.join(outdir, 'beastcmd.log') ]
    return env.Command(target,
                       c['config_beast'],
                       # [ "srun --chdir={} --output=srun.log beast -overwrite -beagle {} >${{TARGETS[3]}} 2>&1".format(outdir, os.path.abspath(str(c['config_beast']))),
                       [ "srun --time=30 --chdir={outdir} --output={outdir}/beastcmd.log beast -overwrite -beagle beast_in.xml >{outdir}/srun.log 2>&1".format(outdir=outdir),
                         Wait(target)
                       ])

@w.add_target_with_env(env)
def mcc(env, outdir, c):
    return env.Command(os.path.join(outdir, 'mcc.tree'),
                        c['runbeast'][2],
                       'treeannotator ${SOURCES} >${TARGET} ')


    

