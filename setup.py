import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="hica",
    version="0.0.2",
    author="Nicolas Sompairac",
    author_email="nicolas.sompairac@curie.fr",
    description="Hierarchical ICA algorithm and applications to genomic data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/NicolasSompairac/Component_Forest",
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
