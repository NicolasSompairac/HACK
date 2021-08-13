import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="hack",
    version="0.0.4",
    author="Nicolas Sompairac",
    author_email="nicolas.sompairac@curie.fr",
    description="Hierachical Analysis of Component linKs (HACK) algorithm for applications to genomic data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/NicolasSompairac/HACK",
    packages=setuptools.find_packages(),
    install_requires = [
        "scipy == 1.6.3",
        "scanpy >= 1.4.5.post2",
        "pandas == 1.2.4",
        "numpy == 1.20.3",
        "networkx == 2.4",
        "matplotlib == 3.4.1",
        "jupyterthemes == 0.20.0",
        "ipywidgets == 7.4.2",
        "requests == 2.22.0",
        "json2html == 1.3.0",
        "mygene==3.2.2"
        ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU GPL License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
