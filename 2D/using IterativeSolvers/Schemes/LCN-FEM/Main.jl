function LCN(Simplices,Nodes,N_of_Simplices,N_of_Nodes,Mesh2Space,SpaceSize,tau,N_t,v,pathway,U0)

Tid=time();





    


#-------------------------Initialize Matrices------------------------------------------------

	diag1 = 0; diag2 = 1; diag3 = Mesh-2; diag4 = Mesh-1;
	vec1 = ones(SpaceSize)*eps(1/10000);  vec2=vec1[1:SpaceSize-diag2]; vec3 = vec1[1:SpaceSize-diag3]; vec4 = vec1[1:SpaceSize-diag4];

    M=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;
    
    A = spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(-diag2=>vec2)+spdiagm(-diag4=>vec4)+spdiagm(diag4=>vec4);

    M_β=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;

    M_v=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;


    M_L=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;


#--------------------------------------------------------------------------------------------




func(x) =1.0;

M = AssembleMatrix_M(Simplices,Nodes,Mesh2Space, func,M,N_of_Simplices,Quad);
M_β = AssembleMatrix_M_formatFEM(Simplices,Nodes,Mesh2Space,  U0,M_β, N_of_Simplices,Quad ); 
M_v = AssembleMatrix_M(Simplices,Nodes,Mesh2Space, v,M_v, N_of_Simplices,Quad);
A = AssembleMatrix_A(Simplices,Nodes,Mesh2Space,A, N_of_Simplices);
#M_L = AssembleMatrix_Lz(Simplices,Nodes,Mesh2Space,M_L, N_of_Simplices,Quad);


Energy = zeros(21,2);

Energy[1,2] = real( dot(U0,(ϵ*A*U0+M_v*U0+β*M_β*U0/2)) ) 



#------------RHS-----------------------------
b1 = zeros(Complex{Float64},SpaceSize)
b1 = M*U0;

b_tot  =  2im*b1;

#Matris = M*2im+tau*(-ϵ*A-β*M_β-M_v+Omega*M_L);
Matris = M*2im+tau*(-ϵ*A-β*M_β-M_v);

Uhat = Matris\b_tot;


U,Energy = TimeStep(Simplices,Nodes,Mesh2Space,N_of_Simplices,N_of_Nodes,β,tau,A,M_β,M,M_v,Uhat,U0,Energy,ϵ,Quad,v,M_L,pathway);
Tid = time()-Tid;

save(pathway*"LCN_beta_"*string(β)*"_Nt_"*string(N_t)*"_Mesh_"*string(Mesh)*".jld","U",U,"Energy",Energy,"Mesh2Space",Mesh2Space,"beta",β,"Tid",Tid)

    
end
