function IM(Simplices::Array{Int64,2}, 
Nodes::Array{Float64,1},  N_of_Simplices::Int64,
N_of_Nodes::Int64, 
Mesh2Space::Array{Int64,1},
SpaceSize::Int64, 
Quad::Dict{String,Array{Float64,1}},  
tau::Float64, 
N_t::Int64,  
v::Function, 
pathway::String, 
U0::Array{Complex{Float64},1})

Tid = time();

    





# ----------------- INITILIZE MATRICES -----------------


diag1 = 0; diag2 = 1;
vec1 = ones(SpaceSize)*eps(1/10000);  vec2=vec1[1:SpaceSize-diag2];
    
M=spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);

A = spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);

M_v=spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);

M_b = spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);

# --------------------------------------------------------------------





#----------------- ASSEMBLE MATRICES -----------------
func(x) = 1.0;

M = AssembleMatrix_M(Simplices,Nodes,Mesh2Space, func,M,N_of_Simplices,Quad);
M_v = AssembleMatrix_M(Simplices,Nodes,Mesh2Space, v,M_v, N_of_Simplices,Quad);
A=ϵ*AssembleMatrix_A(Simplices,Nodes,Mesh2Space,A, N_of_Simplices);
M_b = AssembleMatrix_M_formatFEM(Simplices,Nodes,Mesh2Space,U0, M_b, N_of_Simplices,Quad)

#---------------------------------------------------


#----------------------Do 200 Energy Measurements----------------------"
N_E = 200;
K = T/tau
if(K>N_E)
	E =[1+K/N_E*i for i =0:N_E]
	E = round.(E)
	Energy = zeros(N_E+1,2);
else
	E = [1+i for i = 0:K]
	E = round.(E)
	Energy = zeros(Int(K)+1,2);
end 

Energy[1,1] = real(ϵ*U0'*(A*U0)+β*U0'*(M_b*U0)/2.0)
Energy[1,2] = 0.0

#------------------------------------------------------------------#



U,Energy = TimeStep_IM_FEM(Simplices,Nodes,Mesh2Space,N_of_Simplices,N_of_Nodes,β,tau,N_t,A,M,M_v,U0,Energy,ϵ,T,θ,SpaceSize,Quad,v,M_b,E);



Tid = time()-Tid;

save(pathway*"IM_beta"*string(β)*"_Nt_"*string(N_t)*"_Mesh_"*string(Mesh)*".jld","U",U,"Energy",Energy,"Mesh2Space",Mesh2Space,"beta",β,"Tid",Tid)

    
    
    
end

