function IM(Simplices,Nodes,N_of_Simplices,N_of_Nodes,Mesh2Space,SpaceSize,tau,N_t,v,pathway,U0)
Tid = time();

# --------- INITILIZE MATRICES -----------------


diag1 = 0; diag2 = 1; diag3 = Mesh-2; diag4 = Mesh-1;
vec1 = ones(SpaceSize)*eps(1/10000000);  vec2=vec1[1:SpaceSize-diag2]; vec3 = vec1[1:SpaceSize-diag3]; vec4 = vec1[1:SpaceSize-diag4];

    M=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;
    
    A = spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(-diag2=>vec2)+spdiagm(-diag4=>vec4)+spdiagm(diag4=>vec4);

    M_β=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;

    M_v=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;


#-------------------ASSEMBLE constant MATRICES-------------------------
func(x) = 1.0;
M = AssembleMatrix_M(Simplices,Nodes,Mesh2Space, func,M,N_of_Simplices,Quad);

M_v = AssembleMatrix_M(Simplices,Nodes,Mesh2Space, v,M_v, N_of_Simplices,Quad);


A=ϵ*AssembleMatrix_A(Simplices,Nodes,Mesh2Space,A, N_of_Simplices);
#------------------------------------------------------------------------#


#------------------------- 20 Energy measurements -----------------------#
Energy = zeros(21,2);
#------------------------------------------------------------------------#

U,Energy = TimeStep_IM(Simplices,Nodes,Mesh2Space,N_of_Simplices,N_of_Nodes,β,N_t,tau,A,M,M_v,U0,Energy,ϵ,T,θ,SpaceSize,Quad,Mesh);


Tid = time()-Tid;



save(pathway*"IM_beta_"*string(β)*"_Nt_"*string(N_t)*"_Mesh_"*string(Mesh)*".jld","U",U,"Energy",Energy,"Mesh2Space",Mesh2Space,"beta",β,"Tid",Tid)
     
end

