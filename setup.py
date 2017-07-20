# from distutils.core import setup
from setuptools import setup

setup(
    name = 'PaSDqc',
    description = "Quality control for single cell whole genome sequencing",
    version = '1.0.0',
    packages = ['PaSDqc'],
    scripts = ['scripts/PaSDqc'],
    requires = [
        'numpy',
        'scipy',
        'matplotlib',
        'seaborn',
        'plotly',
        'astropy',
        'statsmodels',
    ],
    author = "Maxwell A. Sherman",
    author_email = "maxwell_sherman@hms.harvard.edu",
    url = "https://github.com/parklab/PaSD-qc",
    license = 'MIT',
    include_package_data = True,
    package_data = {
        'PaSDqc': ['db/*'],
    }
    # include_package_data=True,
)
