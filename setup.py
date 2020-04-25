import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="biohsmm", 
    version="0.1",
    author="Anders Rasmussen",
    author_email="arasmussen@flatironinstitute.org",
    description="Package for making HsMMs for applications in genomics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/anders-w-rasmussen/biohsmm",
    packages=['biohsmm'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
