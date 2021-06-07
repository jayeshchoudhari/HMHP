from setuptools import setup, find_packages

setup(
        name="syntheticDataGeneration",
        version="0.1",
        author="Jayesh Choudhari",
        author_email="choudhari.jayesh@alumni.iitgn.ac.in",
        description="Synthetic Data Generation Pacakge",
        packages=['syntheticDataGeneration'],
        include_package_data=False,
        install_requires=[pkg.strip('\n') for pkg in open("requirements.txt")]
    )
