import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="conspospy",
    version="0.1.0",
    author="Division of Evolutionary Genetics in NIG",
    author_email="hyamashita@nig.ac.jp",
    description="Python version of conspos package.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yamasampo/conspospy",
    packages=setuptools.find_packages(),
    install_requires=[
        'numpy'
    ],
)
