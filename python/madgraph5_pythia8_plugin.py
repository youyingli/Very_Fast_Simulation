#!/usr/bin/env python

import os

batch_script="""#!/bin/bash

#PATH env setup
VFS_PACKAGE_PATH={0}
HEPTOOL=$VFS_PACKAGE_PATH/HepMCTool

source $VFS_PACKAGE_PATH/env_setup.sh
cd $HEPTOOL/MG5_aMC_v2_6_1/
./bin/mg5_aMC -f {1}
"""

def produce_batch_script (madgraph_script):
    bjob_script = batch_script.format(
            os.environ.get('VFSIM_PACKAGE_PATH'),
            madgraph_script
            )
    return bjob_script

if __name__ == '__main__':
    produce_batch_script('<yourscript>')
