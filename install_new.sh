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
STDLIBSETUP="source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-centos7-gcc12-opt/setup.sh"
PYTHONSETUP="source /cvmfs/sft.cern.ch/lcg/releases/LCG_105/Python/3.9.12/x86_64-centos7-gcc12-opt/Python-env.sh"
echo "#!/bin/bash" > $PACKAGE_PATH/env_setup.sh
chmod 755 $PACKAGE_PATH/env_setup.sh
$STDLIBSETUP && echo $STDLIBSETUP >> $PACKAGE_PATH/env_setup.sh
$PYTHONSETUP && echo $PYTHONSETUP >> $PACKAGE_PATH/env_setup.sh
echo 'export VFSIM_PACKAGE_PATH=$PWD' >> $PACKAGE_PATH/env_setup.sh
echo "[INFO] : Finish environment variables setting."

echo "[INFO] : Start to check the packages needed ..."
echo
cd $BACKUP_PATH

########################M A D G R A P H___P Y T H I A 8___I N T E R F A C E###################

if [ ! -f $INSTALL_PATH/MG5_aMC_v3_5_3/HEPTools/MG5aMC_PY8_interface/MG5aMC_PY8_interface ]; then
    echo "[INFO] : M A D G R A P H___P Y T H I A 8___I N T E R F A C E can not be found."
    echo "[INFO] : Start to download M A D G R A P H___P Y T H I A 8___I N T E R F A C E ..."
    echo

    cd $INSTALL_PATH
    wget https://launchpad.net/mg5amcnlo/3.0/3.5.x/+download/MG5_aMC_v3.5.3.tar.gz
    tar -zxf MG5_aMC_v3.5.3.tar.gz
    rm MG5_aMC_v3.5.3.tar.gz
    cd $INSTALL_PATH/MG5_aMC_v3_5_3
    
    OLD_LHAPDF="# lhapdf_py3 = lhapdf-config"
    NEW_LHAPDF=" lhapdf_py3 = /cvmfs/sft.cern.ch/lcg/releases/MCGenerators/lhapdf/6.5.3-3fa11/x86_64-centos7-gcc12-opt/bin/lhapdf-config"
    sed -e "s@$OLD_LHAPDF@$NEW_LHAPDF@g" $INSTALL_PATH/MG5_aMC_v3_5_3/input/.mg5_configuration_default.txt > \
        $INSTALL_PATH/MG5_aMC_v3_5_3/input/mg5_configuration.txt

    echo
    echo "# Package for associated with Madgraph5" > $INSTALL_PATH/install_list
    echo "install pythia8" >> $INSTALL_PATH/install_list
#    echo "install collier" >> $INSTALL_PATH/install_list
#    echo "install ninja" >> $INSTALL_PATH/install_list
#    echo "# aMCatNLO Tool test" >> $INSTALL_PATH/install_list
#    echo "generate p p > z z [QCD]" >> $INSTALL_PATH/install_list
#    echo "output NLOTEST -nojpeg" >> $INSTALL_PATH/install_list
    echo "exit" >> $INSTALL_PATH/install_list
 
    echo "[INFO] : Start to build M A D G R A P H___P Y T H I A 8___I N T E R F A C E ..."
    echo

    ./bin/mg5_aMC -f $INSTALL_PATH/install_list 2>&1 >/dev/null
    rm -rf NLOTEST

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
if [ ! -f ${INSTALL_PATH}/DelphesHepMC2 ]; then

    echo "[INFO] : D E L P H E S 3 package can not be found."
    echo "[INFO] : Start making soft link for D E L P H E S 3 from cvmfs ..."
    echo

    cd $INSTALL_PATH
    ln -s /cvmfs/sft.cern.ch/lcg/releases/delphes/3.5.0-e0c34/x86_64-centos7-gcc12-opt/bin/DelphesHepMC2 DelphesHepMC2

    echo
    echo "[INFO] : D E L P H E S 3 soft link has been made ..."
    echo
    echo

else
    echo "[INFO] : D E L P H E S 3 soft link has been exist."
    echo
    echo
fi


cd $PACKAGE_PATH
echo "[INFO] : All packages have beed set up. Please source ****env_setup.sh**** !!"
echo
