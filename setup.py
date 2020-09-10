import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scGCN",
    version="0.0.1",
    author="QSong",
    author_email="wasqqdyx@gmail.com",
    description="single cell graph convolutional network",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/QSong-github/scGCN",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
