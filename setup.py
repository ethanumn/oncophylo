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
    packages=find_packages(),
    include_package_data=True,
    cmdclass={
        'install': InstallGitHubDependencies,
    }
)