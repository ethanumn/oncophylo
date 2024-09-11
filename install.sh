#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# cd into install script location
cd $SCRIPT_DIR

# Install the package
pip install -e .

# Run the custom install script
bash $SCRIPT_DIR/scripts/install_dependencies.sh