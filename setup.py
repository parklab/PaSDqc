# from distutils.core import setup
from setuptools import setup

def check_dependencies():
    install_requires = []

    # Assuming the python standard library is installed...
    try:
        import pathlib
    except ImportError:
        install_requires.append('pathlib')

    try:
        import numpy
    except ImportError:
        install_requires.append('numpy')

    try:
        import scipy
    except ImportError:
        install_requires.append('scipy')

    try:
        import matplotlib
    except ImportError:
        install_requires.append('matplotlib')

    try:
        import seaborn
    except ImportError:
        install_requires.append('seaborn')

    try:
        import pandas
    except ImportError:
        install_requires.append('pandas')

    try:
        import statsmodels
    except ImportError:
        install_requires.append('statsmodels')

    try:
        import astropy
    except ImportError:
        install_requires.append('astropy')

    try:
        import plotly
    except ImportError:
        install_requires.append('plotly')

    return install_requires

if __name__ == "__main__":
    install_requires = check_dependencies()

    setup(
        name = 'PaSDqc',
        description = "Quality control for single cell whole genome sequencing",
        version = '1.1.0',
        packages = ['PaSDqc'],
        scripts = ['scripts/PaSDqc'],
        install_requires = install_requires,
        author = "Maxwell A. Sherman",
        author_email = "maxwell_sherman@hms.harvard.edu",
        url = "https://github.com/parklab/PaSDqc",
        license = 'MIT',
        include_package_data = True,
        package_data = {
            'PaSDqc': ['db/*'],
        }
        # include_package_data=True,
    )
