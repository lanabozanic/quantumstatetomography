# Quantum State Tomography
This is a simple library is dedicated to performing tomography on an arbitrary amount of qubits / qudits. It uses the Maximum Likelihood Technique as described in [J.B. Altepeter et al.](http://research.physics.illinois.edu/QI/Photonics/tomography-files/tomo_chapter_2004.pdf) to approximate the density matrix of a given quantum state. It is able to perform QST on an arbitary amount of qubits and qudits. Please note that although it can perform these computations the approximation of quantum states can get worse as the amount of qubits / qudits scales.


# Usage

This library is installable with pip. To import the library, run the following command:

```
pip install quantumstatetomography
```

In the library, there are two main classes that are used to perform quantum tomography. the QubitTomo() class and the QuditTomo() class. These are initialized in the following way:

```
import quantumstatetomography as qst
qubit_obj = qst.QubitTomo(n=2)
```
or
```
qudit_obj = qst.QuditTomo(n=1, dim=3)
```


The class is filled with many functions and attributes that allow you to perform tomography on your data and analyse its results both quantatively and visually. For a more comprehensive tutorial on how to use this library, please check out tutorial.ipynb!

# Importing data

This library supports importing quantum tomography data. It is currently only available for the QubitTomo() class. You can import your data as an excel (.xlsx) file whose columns are the measurements taken and their respective counts. For more information on how to format your data, please refer to example.xlsx!
