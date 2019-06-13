include("../Dependencies.jl")


include("Update_Quad.jl")
include("Start_Quad.jl")
include("../Define_Space.jl")
include("../Define_Space_Neuman.jl")


include("../AssembleMatrix_M.jl")
include("../Assemble_M.jl")
include("../AssembleMatrix_A.jl")

include("../Assemble_Matrix_M_formatFEM.jl")
include("../E_Norm.jl")

include("../phi.jl")
include("../Quadrature_7.jl")
include("../L2_Projection.jl")

include("../Reinit_M.jl")

include("Main_NoRot.jl")


