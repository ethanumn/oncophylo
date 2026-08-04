# OncoPhylo

# OncoPhylo

**OncoPhylo** is a lightweight Python package providing utilities for cancer phylogenetics and phylogenetic analysis. It includes a collection of simple tools for working with phylogenetic trees, evolutionary relationships, and related data structures commonly encountered in computational oncology research.

## Features

* Read and manipulate phylogenetic trees
* Compute common tree statistics and metrics
* Utilities for tree traversal and visualization
* Helper functions for evolutionary and cancer genomics analyses
* Lightweight and easy to integrate into existing analysis pipelines

## Installation

Clone the repository:

```bash
git clone https://github.com/ethanumn/oncophylo.git
cd oncophylo
```

Install the package:

```bash
pip install -e .
```

or

```bash
pip install .
```

To automatically install additional tools, see the scripts/ directory. Any tool that depends on an external binary expects the executable to be placed in the oncophylo/bin directory. If the directory does not already exist, you can create it with:

```
mkdir -p oncophylo/bin
```

Ensure that the bin folder is at the same directory level as /ul, /tl, etc.

If you are having trouble install a binary, check scripts/install_dependencies.sh to see how tools would automatically be installed.

## Intended Use

OncoPhylo is designed to provide reusable building blocks for phylogenetic analyses, particularly in cancer genomics research. Rather than being a standalone analysis pipeline, it serves as a collection of utility functions that can be incorporated into larger workflows.

## Requirements

* Python 3.9+
* NumPy
* NetworkX

Additional dependencies are listed in `requirements.txt` or `pyproject.toml`.

