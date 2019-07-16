# Instructions

The code was built for Julia Versions later than 1.0 .

The following packages are needed:

* [IterativeSolvers](https://juliamath.github.io/IterativeSolvers.jl/latest/)
* [SparseArrays](https://docs.julialang.org/en/v1/stdlib/SparseArrays/index.html)
* [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/)
* [JLD](https://github.com/JuliaIO/JLD.jl)
* [FFTW](https://pkg.julialang.org/docs/FFTW)			(For spectral method)

## EXAMPLE 

Start Julia and change directory to Numerical_Comparison/GPE/2D/using IterativeSolvers/Optical_Lattice/, 
this is section 4.1 in the paper. 

Write include("Run.jl"). The files will be saved in the folder T1. 
Once a reference solution has been calculated, the error may be estimated using 
the code in "Error_Est_FEM.jl"



