function RE(Simplices, Nodes,  N_of_Simplices, N_of_Nodes, Mesh2Space, SpaceSize, Quad,  tau, N_t,  v, pathway, U0)



Tid = time();
    
ϕ = zeros(3,N_of_Simplices);


ϕ=Start_Quad(Simplices,Nodes,Mesh2Space,N_of_Simplices,U0,ϕ);
     
    
it = 1;
    
# ---------------------- Initialize Matrices ----------------------

diag1 = 0; diag2 = 1;
vec1 = ones(SpaceSize)*eps(1/10000);  vec2=vec1[1:SpaceSize-diag2];
    
M=spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);
M_ϕ = spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);


func(x) =1.0;
  
M = AssembleMatrix_M(Simplices,Nodes,Mesh2Space, func,M,N_of_Simplices,Quad);

#---------------------- Assemble Matrices ----------------------
A = spzeros(SpaceSize,SpaceSize);
A=spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);
A=AssembleMatrix_A(Simplices,Nodes,Mesh2Space,A, N_of_Simplices);
    
M_v =spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);

M_v = AssembleMatrix_M(Simplices,Nodes,Mesh2Space, v,M_v, N_of_Simplices,Quad);

  

    I = spdiagm(0=>ones(SpaceSize));

  


    U = U0;

    M_b = spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);

#----------------------------------------------------------------------
	

#----------------------Measure energy 200 times----------------------
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

idx_E = 1
Mass = zeros(N_E);
#-----------------------------------------------------------------------------------------------------


for i_t = 1:N_t

	
	t = tau*i_t

	if( i_t == E[idx_E])
		M_b=Reinit_M(Simplices,Nodes,Mesh2Space, M_b, N_of_Simplices);
		M_b = AssembleMatrix_M_formatFEM( Simplices,Nodes,Mesh2Space,U, M_b, N_of_Simplices,Quad)
		Energy[idx_E,2] = real(ϵ*U'*(A*U)+β*U'*(M_b*U)/2.0+U'*M_v*U)
		Energy[idx_E,1] = tau*(i_t-1)
		idx_E +=1;
       end
  
   
        
	ϕ=Update_Quad(Simplices,Nodes,Mesh2Space,N_of_Simplices,U,ϕ);   #Implements rho

        M_ϕ=Reinit_M(Simplices,Nodes,Mesh2Space, M_ϕ, N_of_Simplices);

        M_ϕ = AssembleMatrix_M_RE(Simplices,Nodes,Mesh2Space,ϕ, M_ϕ, N_of_Simplices,Quad);
      
	Matris = M+1im*tau/2*(ϵ*A+β*M_ϕ +M_v)
	RHS = M*U-1im*tau/2*(ϵ*A*U+β*M_ϕ *U+M_v*U);

	U =Tridiagonal(Matris)\RHS;
  


        
      
        
end


#Measure Final Energy
M_b=Reinit_M(Simplices,Nodes,Mesh2Space, M_b, N_of_Simplices);
M_b = AssembleMatrix_M_formatFEM( Simplices,Nodes,Mesh2Space,U, M_b, N_of_Simplices,Quad)
Energy[idx_E,2] = real(ϵ*U'*(A*U)+β*U'*(M_b*U)/2.0+U'*M_v*U)
Energy[idx_E,1] = T

    
Tid = time()-Tid;
save(pathway*"RE_beta"*string(β)*"_Nt_"*string(N_t)*"_Mesh_"*string(N_of_Nodes-1)*".jld","U",U,"Energy",Energy,"Mesh2Space",Mesh2Space,"beta",β,"Tid",Tid)
 
    
    
    
    
    
end






