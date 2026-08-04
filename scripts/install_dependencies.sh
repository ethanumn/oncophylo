#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

LIB_DIR=$SCRIPT_DIR/../oncophylo/lib
BIN_DIR=$SCRIPT_DIR/../oncophylo/bin

# clean
rm -rf $LIB_DIR 
rm -rf $BIN_DIR

# remake directories
mkdir $LIB_DIR 
mkdir $BIN_DIR

# clone all packages from github
echo "Cloning SCITE"
git clone https://github.com/cbg-ethz/SCITE.git $LIB_DIR/SCITE

echo "Cloning infSCITE"
git clone https://github.com/cbg-ethz/infSCITE.git $LIB_DIR/infSCITE

echo "Cloning ConDoR"
git clone https://github.com/raphael-group/ConDoR.git $LIB_DIR/ConDoR
echo "INFO: to use ConDoR you'll need to have GurobiPy properly installed, 'pip install gurobipy', and have Gurobi setup on your system - https://www.gurobi.com/downloads/"

echo "Cloning SPhyR"
git clone https://github.com/elkebir-group/SPhyR.git $LIB_DIR/SPhyR
echo "INFO: to use SPhyR you'll need to install CPLEX (https://www.ibm.com/products/ilog-cplex-optimization-studio)."
echo "Once this is done, you can run the oncophylo/scripts/install_SPhyR.sh script to install LEMON and SPhyR."

echo "Cloning SiFit"
git clone https://github.com/KChen-lab/SiFit.git $LIB_DIR/SiFit
mv $LIB_DIR/SiFit/SiFit.jar $BIN_DIR/SiFit.jar

echo "Cloning HUNTRESS"
git clone https://github.com/PASSIONLab/HUNTRESS.git $LIB_DIR/HUNTRESS

echo "Cloning LoPhy"
git clone https://github.com/ethanumn/LoPhy.git $LIB_DIR/LoPhy
cd LoPhy
mkdir bin
make

# installation for MacOS
if [ "$(uname)" == "Darwin" ]; then

    # download NO-OMP version of COMPASS
    echo "Installing COMPASS"
    cd $LIB_DIR
    curl -LO https://github.com/cbg-ethz/COMPASS/archive/refs/heads/no_OMP.zip
    unzip no_OMP.zip
    cd COMPASS-no_OMP
    make

    # install SCITE
    echo "Installing SCITE"
    cd $LIB_DIR/SCITE
    clang++ *.cpp -o $BIN_DIR/SCITE   

    # install infSCITE
    echo "Installing infSCITE"
    cd $LIB_DIR/infSCITE
    clang++ *.cpp -o $BIN_DIR/infSCITE 


elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then # installation for Linux

    # install linux compatible version of COMPASS
    echo "Cloning COMPASS"
    git clone https://github.com/cbg-ethz/COMPASS $LIB_DIR/COMPASS
    cd $LIB_DIR/COMPASS
    make

    # install SCITE
    echo "Installing SCITE"
    cd $LIB_DIR/SCITE
    g++ *.cpp -o $BIN_DIR/SCITE

    # install infSCITE
    echo "Installing infSCITE"
    cd $LIB_DIR/infSCITE
    g++ *.cpp -o $BIN_DIR/infSCITE

    # install COMPASS
    echo "Installing COMPASS"


else 
    echo "Unknown machine architecture, unable to install dependencies!"
fi