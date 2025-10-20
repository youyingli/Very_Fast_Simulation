#!/usr/bin/env python

from optparse import OptionParser
import os, sys
import random


run_script="""#!/bin/bash

#PATH env setup
VFS_PACKAGE_PATH={0}
HEPTOOL=$VFS_PACKAGE_PATH/HepMCTool/MG5_aMC_v2_9_21/HEPTools

source $VFS_PACKAGE_PATH/env_setup.sh
PYTHIA8=$HEPTOOL/pythia8
PYTHIA8DATA=$HEPTOOL/pythia8/share/Pythia8/xmldoc
SEED=$(($1*$2))

cd $VFS_PACKAGE_PATH/submission/{1}

mkdir -p $1/job-$2
cd $1/job-$2

# Madgraph5
tar zxf {2}
./run.sh {3} $(($1*$2))

#Pythia8
mkdir -p $TEMP/submission/{1}/$1/job-$2
sed -e 's/NEVENTS/{3}/g' {4} > pythia8.dat
sed -i "s|PATH|$TEMP/submission/{1}/$1/job-$2|g" pythia8.dat
sed -i "s/SEED/$SEED/g" pythia8.dat
$HEPTOOL/MG5aMC_PY8_interface/MG5aMC_PY8_interface pythia8.dat

#Delphes
cp {5} .
cp {6} .
sed -i "s/SEED/$SEED/g" {7}
$VFS_PACKAGE_PATH/HepMCTool/DelphesHepMC2 {7} $TEMP/submission/{1}/$1/job-$2/output.root $TEMP/submission/{1}/$1/job-$2/tag_1_pythia8_events.hepmc

#Extraction of jet feature
#root -l -q "$VFS_PACKAGE_PATH/script/event_Extractor.C(\\"$TEMP/submission/{1}/$1/job-$2/output.root\\")"
root -l -q "$VFS_PACKAGE_PATH/script/jetExtractor_modify.C(\\"$TEMP/submission/{1}/$1/job-$2/output.root\\")"
mv output_jet.root $TEMP/submission/{1}/$1/output-$2.root

# Save storgae
cd ../
rm -rf job-$2
rm -rf $TEMP/submission/{1}/$1/job-$2
#rm events.lhe.gz
#rm tag_1_pythia8_events.hepmc

"""

run_script_nlo="""#!/bin/bash

#PATH env setup
VFS_PACKAGE_PATH={0}
HEPTOOL=$VFS_PACKAGE_PATH/HepMCTool/MG5_aMC_v2_9_21/HEPTools

source $VFS_PACKAGE_PATH/env_setup.sh
PYTHIA8=$HEPTOOL/pythia8
PYTHIA8DATA=$HEPTOOL/pythia8/share/Pythia8/xmldoc
SEED=$(($1*$2))

cd $VFS_PACKAGE_PATH/submission/{1}

mkdir -p $1/job-$2
cd $1/job-$2

# Madgraph5
tar zxf {2}
sed -i 's/NEVENTS/{3}/g' Cards/run_card.dat
sed -i "s/SEED/$SEED/g" Cards/run_card.dat
./bin/generate_events --nocompile --only_generation
mv Events/run_01/events.lhe.gz .
rm -rf SubProcesses # Clean Up.

mkdir -p $TEMP/submission/{1}/$1/job-$2

## PileUp ##
sed -e "s/SEED/$SEED/g" $VFS_PACKAGE_PATH/Cards/pythia8/TuneCP5_MinBias.dat > pythia8_MinBias.dat
$VFS_PACKAGE_PATH/HepMCTool/MinBiasGen pythia8_MinBias.dat $TEMP/submission/{1}/$1/job-$2/MinBias.hepmc
hepmc2pileup $TEMP/submission/{1}/$1/job-$2/MinBias.pileup $TEMP/submission/{1}/$1/job-$2/MinBias.hepmc
rm $TEMP/submission/{1}/$1/job-$2/MinBias.hepmc

#Pythia8
sed -e 's/NEVENTS/{3}/g' {4} > pythia8.dat
sed -i "s|PATH|$TEMP/submission/{1}/$1/job-$2|g" pythia8.dat
sed -i "s/SEED/$SEED/g" pythia8.dat
$HEPTOOL/MG5aMC_PY8_interface/MG5aMC_PY8_interface pythia8.dat

#Delphes
cp {5} .
cp {6} .
sed -i "s/SEED/$SEED/g" {7}
sed -i "s|PILEUP|$TEMP/submission/{1}/$1/job-$2/MinBias.pileup|g" {7}
$VFS_PACKAGE_PATH/HepMCTool/DelphesHepMC2 {7} $TEMP/submission/{1}/$1/job-$2/output.root $TEMP/submission/{1}/$1/job-$2/tag_1_pythia8_events.hepmc

#Extraction of jet feature
root -l -q "$VFS_PACKAGE_PATH/script/event_Extractor.C(\\"$TEMP/submission/{1}/$1/job-$2/output.root\\")"
mv output_events.root $TEMP/submission/{1}/$1/output-$2.root
#root -l -q "$VFS_PACKAGE_PATH/script/jetExtractor_modify.C(\\"$TEMP/submission/{1}/$1/job-$2/output.root\\")"
#mv output_jet.root $TEMP/submission/{1}/$1/output-$2.root

# Save storgae
cd ../
rm -rf job-$2
rm -rf $TEMP/submission/{1}/$1/job-$2
#rm events.lhe.gz
#rm tag_1_pythia8_events.hepmc

"""

HTCondorConfig="""executable  = {0}/runjobs.sh
arguments   = $(ClusterId) $(ProcId)

output      = {0}/output/runjob.$(ClusterId).$(ProcId).out
error       = {0}/error/runjob.$(ClusterId).$(ProcId).err
log         = {0}/log/htc.log

#request_cpus = 4
#request_memory = 24GB
#+JobFlavour = "large"

max_materialize = 20
max_retries = 1
queue {1}
"""


def Option_Parser(argv):

    usage='usage: %prog [options] arg\n'
    usage+='The script handles the MC generation by parallel calculation.\n'
    usage+='It can automatically produce some scripts for given MCTool and submit to batch system.\n'
    usage+='For more information, please see README !!!\n'
    parser = OptionParser(usage=usage)

    parser.add_option('-t', '--tag',
            type='str', dest='tag', default='Standard',
            help='tag output MC samples'
            )
    parser.add_option('--madgraph_gridpack',
            type='str', dest='madgraph_gridpack', default='',
            help='madgraph gridpack path'
            )
    parser.add_option('-n', '--nEvents',
            type='int', dest='nEvents', default=10000,
            help='Number of events per job'
            )
    parser.add_option('-j', '--nJobs',
            type='int', dest='nJobs', default=10,
            help='Number of jobs'
            )
    parser.add_option('--pythia8_card',
            type='str', dest='pythia8_card', default=str(os.environ.get('VFSIM_PACKAGE_PATH')) + '/Cards/pythia8/TuneCP5.dat',
            help='Datacard needed when running the pythia8'
            )
    parser.add_option('--delphes_card',
            #type='str', dest='delphes_card', default=str(os.environ.get('VFSIM_PACKAGE_PATH')) + '/Cards/Delphes/delphes_card_CMS.tcl',
            #type='str', dest='delphes_card', default=str(os.environ.get('VFSIM_PACKAGE_PATH')) + '/Cards/Delphes/delphes_card_CMS_PileUp_Puppi.tcl',
            type='str', dest='delphes_card', default=str(os.environ.get('VFSIM_PACKAGE_PATH')) + '/Cards/Delphes/delphes_card_CMS_PhaseII_Puppi.tcl',
            help='Datacard needed when running the detector simulation by Delphes3'
            )
    parser.add_option('--delphes_card_support',
            #type='str', dest='delphes_card_support', default=str(os.environ.get('VFSIM_PACKAGE_PATH')) + '/Cards/Delphes/trackResolutionCMS.tcl',
            type='str', dest='delphes_card_support', default=str(os.environ.get('VFSIM_PACKAGE_PATH')) + '/Cards/Delphes/muonMomentumResolution.tcl',
            help='Supported datacard needed when running the detector simulation by Delphes3'
            )
    parser.add_option('--isNLO',
            action='store_true', dest='isNLO',
            help='Is NLO'
            )

    (options, args) = parser.parse_args(argv)
    return options


def HTCondor (argv):

    options = Option_Parser(argv)

    package_path = os.environ.get('VFSIM_PACKAGE_PATH')
    job_path = f'{package_path}/submission/{options.tag}'
    os.system(f'mkdir -p {job_path}')


    with open(f'{job_path}/runjobs.sh', 'w') as job_script :

        job_script.write(globals()['run_script_nlo' if options.isNLO else 'run_script'].format(
                            package_path,
                            options.tag,
                            options.madgraph_gridpack,
                            options.nEvents,
                            options.pythia8_card,
                            options.delphes_card,
                            options.delphes_card_support,
                            options.delphes_card[options.delphes_card.rfind('/')+1:]
                                            ))
        job_script.write('\n')


    with open(f'{job_path}/runjobs.sub', 'w') as condor_config :

        condor_config.write(HTCondorConfig.format(
                            job_path,
                            options.nJobs
                                            ))
        condor_config.write('\n')

        os.system(f'mkdir -p {job_path}/error')
        os.system(f'mkdir -p {job_path}/output')
        os.system(f'mkdir -p {job_path}/log')

        print( f'{job_path}/runjobs.sub' )

    os.system(f'condor_submit {job_path}/runjobs.sub')


if  __name__ == '__main__':
    sys.exit( HTCondor(sys.argv[1:]) )
