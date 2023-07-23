import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="hack",
    version="0.0.5",
    author="Nicolas Sompairac",
    author_email="nicolas.sompairac@curie.fr",
    description="Hierachical Analysis of Component linKs (HACK) algorithm for applications to genomic data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/NicolasSompairac/HACK",
    packages=setuptools.find_packages(),
    install_requires = [
        "scipy == 1.9.1",
        "scanpy >= 1.9.1",
        "stabilized-ica == 2.0.0",
        "pandas == 1.5.0",
        "numpy == 1.23.3",
        "networkx == 2.8.6",
        "matplotlib == 3.6.0",
        "jupyterthemes == 0.20.0",
        "ipywidgets == 7.4.2",
        "requests == 2.28.1",
        "json2html == 1.3.0",
        "mygene==3.2.2"
        ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU GPL License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8.13',
)
