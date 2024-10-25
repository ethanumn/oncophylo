import os, subprocess
from setuptools import setup, find_packages
from setuptools.command.install import install as setuptools_install

class InstallGitHubDependencies(setuptools_install):
    """Install GitHub dependencies when installing OncoPhylo"""
    
    def run(self):
        # Run the original install command
        setuptools_install.run(self)
        
        # Now run the custom dependencies installation script
        script_path = os.path.join(os.path.dirname(__file__), 'install_dependencies.sh')
        if os.path.exists(script_path):
            subprocess.check_call(['bash', script_path])
        else:
            print(f"Warning: {script_path} does not exist.")
        
setup(
    name='OncoPhylo',
    version='0.1',
    description = "A package to perform cancer phylogeny reconstruction.",
    python_requires = ">=3.8",
    packages=find_packages(),
    include_package_data=True,
    cmdclass={
        'install': InstallGitHubDependencies,
    },
    install_requires=[
            "setuptools>=61.0", 
            "numpy>=1.21.0",
            "pandas>=2.0",
            "seaborn",
            "scipy>=1.1.0", 
            "kmodes>=0.12", 
            "networkx>=3.0",
            "pydot>=2.0",
            "gurobipy",
            "anndata"
    ]
)