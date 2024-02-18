#!/usr/bin/env python

import os

batch_script="""#!/bin/bash

#PATH env setup
VFS_PACKAGE_PATH={0}
HEPTOOL=$VFS_PACKAGE_PATH/HepMCTool/MG5_aMC_v3_5_3/HEPTools

source $VFS_PACKAGE_PATH/env_setup.sh
PYTHIA8=$HEPTOOL/pythia8
PYTHIA8DATA=$HEPTOOL/pythia8/share/Pythia8/xmldoc


mkdir job-$1
cd job-$1

# Madgraph5
ln -s {1} run.sh
./run.sh {2} $1

#Pythia8
$HEPTOOL/MG5aMC_PY8_interface/MG5aMC_PY8_interface {3}

#Delphes
$VFS_PACKAGE_PATH/HepMCTool/DelphesHepMC2 {4} output-$1.root tag_1_pythia8_events.hepmc

# Save storgae
rm events.lhe.gz
rm tag_1_pythia8_events.hepmc


"""

def produce_batch_script (madgraph_script):
    bjob_script = batch_script.format(
            os.environ.get('VFSIM_PACKAGE_PATH'),
            madgraph_gridpack,
            nEvents,
            pythia8_card,
            delphes_card,
            )
    return bjob_script

if __name__ == '__main__':
    produce_batch_script('<yourscript>')
