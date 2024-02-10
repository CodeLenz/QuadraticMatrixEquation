# QuadraticMatrixEquation

[![Build Status](https://github.com/eduardobarplenz@gmail.com/QuadraticMatrixEquation/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/eduardobarplenz@gmail.com/QuadraticMatrixEquation/actions/workflows/CI.yml?query=branch%3Amain)

Solve 

$AX² + BX + C = 0$

and

$X²A + XB + C = 0$

where $A$, $B$ and $C$ are $n \times n$ complex matrices.

using Newton-Raphson.

Observation: The Sylvester equation to find the Fréchet derivative is evaluated using the MatrixEquations package (https://github.com/andreasvarga/MatrixEquations.jl)
