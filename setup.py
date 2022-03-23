import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="qubit_tomography-lanabozanic",
    version="0.0.1",
    author="Lana Bozanic",
    author_email="lbozanic@uottawa.ca",
    description="A python package used to perform quantum state tomography.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lanabozanic/qubit_tomography",
    project_urls={
        "Bug Tracker": "https://github.com/lanabozanic/qubit_tomography/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)