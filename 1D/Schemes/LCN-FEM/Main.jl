function LCN(Simplices::Array{Int64,2},
Nodes::Array{Float64,1},
N_of_Simplices::Int64,
N_of_Nodes::Int64,
Mesh2Space::Array{Int64,1},
SpaceSize::Int64,
Quad::Dict{String,Array{Float64,1}},
tau::Float64,
N_t,
v::Function,
pathway::String,
U0::Array{Complex{Float64},1},)

Tid=time();




    

#------------Initialize Matrices------------


diag1 = 0; diag2 = 1;
vec1 = ones(SpaceSize)*eps(1/10000);  vec2=vec1[1:SpaceSize-diag2];

M=spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);
M_β = spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);
M_v = spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);
A = spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);
#-------------------------------------------

#------------ Assemble Matrices ------------
func(x) =1.0;

M = AssembleMatrix_M(Simplices,Nodes,Mesh2Space, func,M,N_of_Simplices,Quad);
M_β = AssembleMatrix_M_formatFEM(Simplices,Nodes,Mesh2Space,  U0,M_β, N_of_Simplices,Quad ); 
M_v = AssembleMatrix_M(Simplices,Nodes,Mesh2Space, v,M_v, N_of_Simplices,Quad);
A=AssembleMatrix_A(Simplices,Nodes,Mesh2Space,A, N_of_Simplices);

#------------------------------------------------



#------------------------RHS------------------------
b = M*U0;
b = 2im*b;
#---------------------------------------------------


Matris = M*2im+tau*(-ϵ*A-β*M_β-M_v);

Uhat = Tridiagonal(Matris)\b;

U,Energy = TimeStep_LCN(Simplices,Nodes,Mesh2Space,N_of_Simplices,N_of_Nodes,β,tau,N_t,A,M_β,M,M_v,Uhat,U0,ϵ,Quad,v);
Tid = time()-Tid;

save(pathway*"LCN_beta"*string(β)*"_Nt_"*string(N_t)*"_Mesh_"*string(N_of_Nodes-1)*".jld","U",U,"Energy",Energy,"Mesh2Space",Mesh2Space,"beta",β,"Tid",Tid)
 
    

end
