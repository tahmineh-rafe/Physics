# Unitary Matrix Decomposition for Photonic Quantum Computing

## Overview

This project implements the decomposition of a unitary matrix into the physical parameters required for its implementation in an integrated photonic circuit.

The decomposition procedure is based on the unitary matrix decomposition methods introduced by Reck et al. and Clements et al. [1, 2, 3]. These methods provide systematic approaches for decomposing an arbitrary unitary matrix into a network of two-mode transformations.

In this project, the decomposition algorithms were implemented independently in Python. To test the implementation, a unitary matrix obtained from the work of Arrazola et al. was used as the target matrix [4].

## What does the code do?

The code takes a target unitary matrix as its input and decomposes it into a set of parameters, including mixing angles (θ) and phase parameters (φ).

These parameters describe the two-mode transformations required to construct the target unitary transformation.

The extracted parameters are then used to determine the arrangement and parameters of Mach–Zehnder interferometers (MZIs) in a photonic circuit.

The resulting circuit is designed such that, when the calculated θ and φ parameters are applied to the corresponding MZIs, the original target unitary matrix can be reconstructed.

In other words, the computational procedure establishes a connection between:

**Target unitary matrix → Decomposition → θ and φ parameters → Photonic circuit → Reconstructed unitary matrix**

## Physical interpretation

A unitary matrix describes a transformation that can be implemented in a lossless linear optical system.

In an integrated photonic circuit, the required unitary transformation can be implemented using elementary optical components such as beam splitters, phase shifters, and Mach–Zehnder interferometers.

The decomposition process therefore provides a way to translate a mathematically specified unitary transformation into experimentally realizable photonic circuit parameters.

The calculated θ and φ parameters determine how the individual Mach–Zehnder interferometers should be configured so that the complete photonic network performs the desired unitary transformation.

## Application in Quantum Computing

This approach is particularly relevant to photonic quantum computing.

Quantum gates and more general quantum operations can be represented mathematically by unitary transformations. Therefore, if a desired quantum operation is represented by a unitary matrix, unitary decomposition methods can be used to determine how that transformation can be implemented using a physical photonic circuit.

In this context, the decomposition provides a bridge between the abstract mathematical description of a quantum operation and its physical implementation on an integrated photonic platform.

The overall concept can be summarized as:

**Quantum operation (unitary matrix)**
↓
**Unitary matrix decomposition**
↓
**Extraction of θ and φ parameters**
↓
**Configuration of Mach–Zehnder interferometers**
↓
**Implementation of the desired transformation in a photonic circuit**

## Purpose of this Repository

The purpose of this repository is to provide a Python implementation of unitary matrix decomposition and to demonstrate how the resulting decomposition parameters can be related to the implementation of the corresponding unitary transformation using Mach–Zehnder interferometer-based photonic networks.

The repository focuses on the computational implementation, matrix reconstruction, and numerical verification of the decomposition procedure.

## Citation

If you use this implementation in academic work, please cite this repository as:

Rafe, Tahmineh. (2026).
Unitary Matrix Decomposition for Photonic Quantum Computing.
GitHub repository: Physics.

## References

[1] Reck, M., Zeilinger, A., Bernstein, H. J., & Bertani, P. (1994). Experimental realization of any discrete unitary operator. Physical review letters, 73(1), 58.

[2] Clements, W. R., Humphreys, P. C., Metcalf, B. J., Kolthammer, W. S., & Walmsley, I. A. (2016). Optimal design for universal multiport interferometers. Optica, 3(12), 1460-1465.

[3] Capmany, J., & Pérez, D. (2020). Programmable integrated photonics. Oxford University Press.

[4] Arrazola, J. M., Bergholm, V., Brádler, K., Bromley, T. R., Collins, M. J., Dhand, I., ... & Zhang, Y. (2021). Quantum circuits with many photons on a programmable nanophotonic chip. Nature, 591(7848), 54-60.
