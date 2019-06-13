using IterativeSolvers
using SparseArrays
using LinearAlgebra
using JLD


#---------------Load Functions in common --------------------------
include("Dirichlet_BC.jl")

include("AssembleMatrix_M_formatFEM.jl")
include("AssembleMatrix_M.jl")
include("AssembleMatrix_A.jl")
include("phi.jl")
include("Reinit_M.jl")
include("Quadrature_7.jl")

include("L2_Projection.jl")

#------------------------------------------------------------------



#---------------------Load All schemes-------------------------------
include("./CN-FEM/CN_v01.jl")
include("./IM-FEM/IM_v01.jl")
include("./LCN-FEM/LCN_v01.jl")
include("./RE-FEM/RE_v01.jl")
#--------------------------------------------------------------------
