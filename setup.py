
from setuptools import setup, find_packages

NAME = 'cabean'

setup(name=NAME,
    version='9999',
    author = "Loïc Paulevé",
    author_email = "loic.pauleve@labri.fr",
    url = "https://github.com/algorecell/pyCabean",
    description = "Python interface to CABEAN, A Software Tool for the Control of Asynchronous Boolean Networks",
    install_requires = [
        "colomoto_jupyter",
        "algorecell_types",
    ],
    classifiers=[
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords="computational systems biology",

    include_package_data = True,
    packages = find_packages(),
    py_modules = ["cabean_setup"]
)

