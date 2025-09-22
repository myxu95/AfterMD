from setuptools import setup, find_packages

setup(
    name="aftermd",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "MDAnalysis>=2.6.0",
        "numpy>=1.24.0",
        "matplotlib>=3.7.0",
        "seaborn>=0.12.0",
        "pandas>=2.0.0",
        "scipy>=1.10.0",
        "plotly>=5.15.0",
    ],
    author="Research Team",
    description="GROMACS MD analysis toolkit",
    python_requires=">=3.8",
)