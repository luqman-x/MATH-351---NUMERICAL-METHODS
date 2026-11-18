# Numerical Methods

This repository contains Python examples and a mathematical derivation for numerical methods used to solve systems of linear equations. The examples focus on Gaussian elimination, pivoting, and LU decomposition, with step-by-step console output in several scripts.

## Contents

| File | Description |
| --- | --- |
| [`Naive_Gaussian_Elimination_Algorithm.py`](Naive_Gaussian_Elimination_Algorithm.py) | Basic Gaussian elimination with forward elimination and back substitution. |
| [`Step_by_Step_Process.py`](Step_by_Step_Process.py) | Naive Gaussian elimination with detailed printed steps, including row operations and intermediate matrices. |
| [`Partial_Pivoting.py`](Partial_Pivoting.py) | Gaussian elimination with partial pivoting, choosing the largest available pivot in each column. |
| [`LU_Decomposition.py`](LU_Decomposition.py) | LU decomposition with partial pivoting, verification of $PA=LU$, and forward and backward substitution to solve the system. |
| [`Jacobi iteration mathematical derivation/README.md`](Jacobi%20iteration%20mathematical%20derivation/README.md) | Derives a scalar recurrence and direct formulas for selected iterations of a two-variable Jacobi method, with a worked example. |

The Python scripts include example systems in the source code and print or return their computed solutions. They are intended as instructional examples; the naive elimination implementation does not use pivoting, while the pivoting examples are more robust to a zero or small pivot.

## Requirements

- Python 3
- NumPy

Install NumPy with:

```bash
python -m pip install numpy
```

## Run an example

From the repository folder, run a script with Python. For example:

```bash
python Partial_Pivoting.py
```

Other scripts can be run the same way:

```bash
python Naive_Gaussian_Elimination_Algorithm.py
python Step_by_Step_Process.py
python LU_Decomposition.py
```

Each script currently runs its included example when executed. To solve a different system, edit the coefficient matrix `A` and right-hand-side vector `B` in that script.

## Methods at a glance

- **Gaussian elimination** transforms the system to upper-triangular form, then uses back substitution.
- **Partial pivoting** swaps rows to use a larger-magnitude pivot, which helps avoid division by zero and can improve numerical stability.
- **LU decomposition** factors a permuted matrix as $PA=LU$, then solves using forward and backward substitution.
- **Jacobi iteration** is covered in the derivation document. That document shows how a two-variable iteration can be expressed as a scalar recurrence and how even or odd indexed values can be computed directly.

## Notes

These programs are small educational implementations rather than a general-purpose linear algebra package. In practice, use established numerical libraries for production calculations, especially for large or poorly conditioned systems.
