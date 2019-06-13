function TWOSTEP(Simplices, Nodes,  N_of_Simplices, N_of_Nodes, Mesh2Space, SpaceSize, Quad,  tau, N_T,  v, pathway, U0)

Tid = time()

    w0 = U0;
    w1 = U0;

#-----------------------------Do 200 Energy Measurments ------------------------------------------"
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
#-----------------------------------------------------------------------------------#

    
    #q = -β;
    #Initialize Matrices

	diag1 = 0; diag2 = 1;
	vec1 = ones(SpaceSize)*eps(1/10000);  vec2=vec1[1:SpaceSize-diag2];
    
 	M=spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);
	A=spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);
	M_v=spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);	
    	M_β = spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);
	I = spdiagm(0=>ones(SpaceSize));
   
    func(x) =1.0;
    M = AssembleMatrix_M(Simplices,Nodes,Mesh2Space, func,M,N_of_Simplices,Quad);

    #----------Create A matrix <du,dv>----------------------------
    A = ϵ*AssembleMatrix_A(Simplices,Nodes,Mesh2Space,A, N_of_Simplices);
    
    M_v = AssembleMatrix_M(Simplices,Nodes,Mesh2Space, v,M_v, N_of_Simplices,Quad);


    
    # -----------------------------------------------------------------

    
 

    #STEP ONE
    M_β = AssembleMatrix_M_formatFEM(Simplices,Nodes,Mesh2Space,  w0,M_β, N_of_Simplices,Quad);
    Matris = (M+1im*tau/4*(A+M_v+β*M_β));
    RHS = M*w0-1im*tau/4*(A*w0+M_v*w0+β*M_β*w0);
    w05 = Matris\RHS;

#Measure initial energy
	Energy[idx_E,2] = real(w0'*(A*w0)+β*w0'*(M_β*w0)/2.0+w0'*M_v*w0)[1]
        Energy[idx_E,1] = 0.

    #STEP TWO
    M_β=Reinit_M(Simplices,Nodes,Mesh2Space, M_β, N_of_Simplices);
    M_β=AssembleMatrix_M_formatFEM(Simplices,Nodes,Mesh2Space,  w05[:],M_β, N_of_Simplices,Quad);
    Matris = (M+1im*tau/2*(A+M_v+β*M_β));
    RHS = M*w0-1im*tau/2*(A*w0+M_v*w0+β*M_β*w0);
    w1 = Matris\RHS;

idx_E = 2;


    for idx_t = 2:N_t

 
      
    M_β=Reinit_M(Simplices,Nodes,Mesh2Space, M_β, N_of_Simplices);
    M_β=AssembleMatrix_M_formatFEM(Simplices,Nodes,Mesh2Space,  w1[:],M_β, N_of_Simplices,Quad);

      if( idx_t == E[idx_E]) 
	    Energy[idx_E,2] = real(w1'*(A*w1)+β*w1'*(M_β*w1)/2.0+w1'*M_v*w1)[1]
	    Energy[idx_E,1] = (idx_t-1)*tau
	    idx_E +=1;	

       end
	

       Matris = M+1im*tau*(A+M_v+M_β*β);
   
        RHS = M*w0-1im*tau*(A*w0+M_v*w0+β*M_β*w0);
        
        w0 = w1;
        w1 = Tridiagonal(Matris)\RHS;
       
	t = tau * idx_t;
        
      
        

    end

M_β=Reinit_M(Simplices,Nodes,Mesh2Space, M_β, N_of_Simplices);
M_β=AssembleMatrix_M_formatFEM(Simplices,Nodes,Mesh2Space,  w1[:],M_β, N_of_Simplices,Quad);

Energy[end,2] = real(w1'*(A*w1)+β*w1'*(M_β*w1)/2.0+w1'*M_v*w1)[1]
Energy[end,1] = T

   
Tid = time()-Tid;

 save(pathway*"TWOSTEP_beta"*string(β)*"_Nt_"*string(N_t)*"_Mesh_"*string(N_of_Nodes-1)*".jld","U",w1,"Energy",Energy,"Mesh2Space",Mesh2Space,"beta",β,"Tid",Tid)
 
    

    
end
