#!/bin/bash

echo "********   Welcome to Very_fast_simulation package   ********"

PACKAGE_PATH=$PWD
BACKUP_PATH=$PACKAGE_PATH/BackUp
INSTALL_PATH=$PACKAGE_PATH/HepMCTool

if [ ! -d $INSTALL_PATH ] || [ ! -d $BACKUP_PATH ]; then
    mkdir $BACKUP_PATH
    mkdir $INSTALL_PATH
fi

echo "[INFO] : Setting environment variables ..."
STDLIBSETUP="source /cvmfs/sft.cern.ch/lcg/releases/LCG_87/gcc/4.9.3/x86_64-slc6/setup.sh"
ROOTLIBSETUP="source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.36/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh"
PYTHONSETUP="source /opt/rh/python27/enable"
echo "#!/bin/bash" > $PACKAGE_PATH/env_setup.sh
chmod 755 $PACKAGE_PATH/env_setup.sh
$STDLIBSETUP && echo $STDLIBSETUP >> $PACKAGE_PATH/env_setup.sh
$ROOTLIBSETUP && echo $ROOTLIBSETUP >> $PACKAGE_PATH/env_setup.sh
$PYTHONSETUP && echo $PYTHONSETUP >> $PACKAGE_PATH/env_setup.sh
echo 'export VFSIM_PACKAGE_PATH=$PWD' >> $PACKAGE_PATH/env_setup.sh
echo "[INFO] : Finish environment variables setting."

echo "[INFO] : Start to check the packages needed ..."
echo
cd $BACKUP_PATH


#########################L H A P D F 6 . 1 . 6#################################################
if [ ! -f $INSTALL_PATH/lhapdf_install/bin/lhapdf-config ]; then

    echo "[INFO] : L H A P D F 6 . 1 . 6 plugin can not be found."
    echo "[INFO] : Start to download L H A P D F 6 . 1 . 6 plugin ..."
    echo

    wget http://www.hepforge.org/archive/lhapdf/LHAPDF-6.1.6.tar.gz
    tar -zxf LHAPDF-6.1.6.tar.gz
    rm LHAPDF-6.1.6.tar.gz
    cd $BACKUP_PATH/LHAPDF-6.1.6

    echo
    echo "[INFO] : Start to build L H A P D F 6 . 1 . 6 plugin ..."
    echo

    ./configure --with-boost=/afs/cern.ch/sw/lcg/external/Boost/1.55.0_python2.7/x86_64-slc6-gcc47-opt \
        --prefix=${INSTALL_PATH}/lhapdf_install 2>&1 >/dev/null
    make -j 10 && make install 2>&1 >/dev/null
    $INSTALL_PATH/lhapdf_install/bin/lhapdf install NNPDF23_lo_as_0130_qed

    echo
    echo "[INFO] : L H A P D F 6 . 1 . 6 plugin has been finished."
    echo
    echo


else

    echo "[INFO] : L H A P D F 6 . 1 . 6 plugin has been exist."
    echo
    echo
fi

echo 'export PATH=$PWD/HepMCTool/lhapdf_install/bin:$PATH' >> $PACKAGE_PATH/env_setup.sh
echo 'export LD_LIBRARY_PATH=$PWD/HepMCTool/lhapdf_install/lib:$LD_LIBRARY_PATH' >> $PACKAGE_PATH/env_setup.sh
echo 'export PYTHONPATH=$PWD/HepMCTool/lhapdf_install/lib64/python2.7/site-packages:$PYTHONPATH' >> $PACKAGE_PATH/env_setup.sh


########################M A D G R A P H___P Y T H I A 8___I N T E R F A C E###################

if [ ! -f $INSTALL_PATH/MG5_aMC_v2_6_1/HEPTools/MG5aMC_PY8_interface/MG5aMC_PY8_interface ]; then
    echo "[INFO] : M A D G R A P H___P Y T H I A 8___I N T E R F A C E can not be found."
    echo "[INFO] : Start to download M A D G R A P H___P Y T H I A 8___I N T E R F A C E ..."
    echo

    cd $INSTALL_PATH
    wget https://launchpad.net/mg5amcnlo/2.0/2.6.x/+download/MG5_aMC_v2.6.1.tar.gz
    tar -zxf MG5_aMC_v2.6.1.tar.gz
    rm MG5_aMC_v2.6.1.tar.gz
    cd $INSTALL_PATH/MG5_aMC_v2_6_1
    
    OLD_LHAPDF="# lhapdf = lhapdf-config"
    NEW_LHAPDF=" lhapdf = ${INSTALL_PATH}/lhapdf_install/bin/lhapdf-config"
    sed -e "s@$OLD_LHAPDF@$NEW_LHAPDF@g" $INSTALL_PATH/MG5_aMC_v2_6_1/input/.mg5_configuration_default.txt > \
        $INSTALL_PATH/MG5_aMC_v2_6_1/input/mg5_configuration.txt
    
    echo
    echo "# Package for associated with Madgraph5" > $INSTALL_PATH/install_list
    echo "install pythia8" >> $INSTALL_PATH/install_list
    echo "exit" >> $INSTALL_PATH/install_list
 
    echo "[INFO] : Start to build M A D G R A P H___P Y T H I A 8___I N T E R F A C E ..."
    echo

    ./bin/mg5_aMC -f $INSTALL_PATH/install_list 2>&1 >/dev/null

    echo
    echo "[INFO] : M A D G R A P H___P Y T H I A 8___I N T E R F A C E package has been finished."
    echo
    echo

else
    echo "[INFO] : M A D G R A P H___P Y T H I A 8___I N T E R F A C E has been exist."
    echo
    echo
fi


##########################D E L P H E S 3#######################################################


if [ ! -f ${INSTALL_PATH}/delphes_install/bin/DelphesHepMC ]; then

    echo "[INFO] : D E L P H E S 3 package can not be found."
    echo "[INFO] : Start to download package D E L P H E S 3 ..."
    echo

    cd $BACKUP_PATH
    git clone https://github.com/delphes/delphes.git
    cd $BACKUP_PATH/delphes
    mkdir $BACKUP_PATH/delphes/build
    cd $BACKUP_PATH/delphes/build
    
    echo
    echo "[INFO] : Start to build D E L P H E S 3 package ..."
    echo
    
    cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_PATH}/delphes_install .. 2>&1 >/dev/null
    make -j 4 install 2>&1 >/dev/null
    
    echo
    echo "[INFO] : D E L P H E S 3 package has been finished..."
    echo
    echo

else
    echo "[INFO] : D E L P H E S 3 package has been exist."
    echo
    echo
fi

cd $PACKAGE_PATH
echo "[INFO] : All packages have beed set up. Please source ****env_setup.sh**** !!"
echo
echo
