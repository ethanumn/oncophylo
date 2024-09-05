#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

LIB_DIR=$SCRIPT_DIR/oncophylo/lib
BIN_DIR=$SCRIPT_DIR/oncophylo/bin

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

if [ "$(uname)" == "Darwin" ]; then

    # install SCITE
    echo "Installing SCITE"
    cd $LIB_DIR/SCITE
    clang++ *.cpp -o $BIN_DIR/SCITE   

    # install infSCITE
    echo "Installing infSCITE"
    cd $LIB_DIR/infSCITE
    clang++ *.cpp -o $BIN_DIR/infSCITE   

elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then

    # install SCITE
    echo "Installing SCITE"
    cd $LIB_DIR/SCITE
    g++ *.cpp -o $BIN_DIR/SCITE

    # install infSCITE
    echo "Installing infSCITE"
    cd $LIB_DIR/infSCITE
    g++ *.cpp -o $BIN_DIR/infSCITE

else 
    echo "Unknown machine architecture, unable to install dependencies!"
fi