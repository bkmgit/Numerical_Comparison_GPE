# Instructions

The code was built for Julia Versions later than 1.0 .

The following packages are needed:


* [SparseArrays](https://docs.julialang.org/en/v1/stdlib/SparseArrays/index.html)
* [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/)
* [JLD](https://github.com/JuliaIO/JLD.jl)
* [FFTW](https://pkg.julialang.org/docs/FFTW)			(For spectral method)

## EXAMPLE 

Start Julia and change directory to Numerical_Comparison/GPE/1D/Single_Soliton, write include("Run.jl"). 
The code solves 
<img src="https://latex.codecogs.com/svg.latex?\Large&space;iu_t=-u_{xx}-|u|^2u" title="\Large Iu_t=-u_{xx}-|u|^2u" />
with initial value 
<img src="https://latex.codecogs.com/svg.latex?\Large&space;u0=\sqrt{2}\textup{exp}\left(\frac{ix}{2}\right)\textup{sech}(x)" title="\Large iu_t=-u_{xx}-|u|^2u" />
which is section 3.1 in the paper. The files are saved in the folder T10. 

To calculate the errors simply run the command include("Error_Est_FEM.jl")




