import os, subprocess
from setuptools import setup, find_packages
from setuptools.command.install import install as setuptools_install

class InstallGitHubDependencies(setuptools_install):
    """Install GitHub dependencies when installing oncophylo"""
    
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
    name='oncophylo',
    version='1.0',
    description = "oncophylo is a package for inferring cancer evolution from sequencing data",
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
            "anndata",
            "rpy2>=3.5.11"
    ]
)