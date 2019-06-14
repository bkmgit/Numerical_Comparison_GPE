function RE(Simplices,Nodes,N_of_Simplices,N_of_Nodes,Mesh2Space,SpaceSize,tau,N_t,v,pathway,U0)

# --------- INITILIZE MATRICES! -----------------

Tid = time();
println(SpaceSize)
	diag1 = 0; diag2 = 1; diag3 = Mesh-2; diag4 = Mesh-1;
	vec1 = ones(SpaceSize)*eps(1/10000);  vec2=vec1[1:SpaceSize-diag2]; vec3 = vec1[1:SpaceSize-diag3]; vec4 = vec1[1:SpaceSize-diag4];

    M=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;
    
A = spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(-diag2=>vec2)+spdiagm(-diag4=>vec4)+spdiagm(diag4=>vec4);


M_ϕ=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;

#M_ϕref = spdiagm(vec1,0)+spdiagm(vec2,diag2,SpaceSize,SpaceSize)+spdiagm(vec3,diag3,SpaceSize,SpaceSize)+spdiagm(vec4,diag4,SpaceSize,SpaceSize)+spdiagm(vec2,-diag2,SpaceSize,SpaceSize)+spdiagm(vec3,-diag3,SpaceSize,SpaceSize)+spdiagm(vec4,-diag4,SpaceSize,SpaceSize) ;

    M_v=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;


    M_β=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;

# -----------------------------------------------

ϕ = zeros(6,N_of_Simplices);
#Phiold = zeros(6,N_of_Simplices);
ϕ=Start_Quad(Simplices,Nodes,Mesh2Space,N_of_Simplices,U0,ϕ);
σ = 1;
      
#-----------------------------Assemble Matrices ------------------------------------------------------------------- 

    func(x) =1.0;

#tic()
    M = AssembleMatrix_M(Simplices,Nodes,Mesh2Space, func,M,N_of_Simplices,Quad);
#toc()

#tic()
    A=AssembleMatrix_A(Simplices,Nodes,Mesh2Space,A, N_of_Simplices);
#toc()

#tic()
    M_v = AssembleMatrix_M(Simplices,Nodes,Mesh2Space,v,M_v, N_of_Simplices,Quad);
#toc()
    M_β = AssembleMatrix_M_formatFEM(Simplices,Nodes,Mesh2Space,  U0[:],M_β, N_of_Simplices,Quad ); 

# -----------------------------------------------------------------

    
    U = U0;
    
#---------------------- Energy ---------------------------------------
Energy = zeros(21,2);
Energy[1,2] = real( dot(U,(0.5*A*U+M_v*U+β*M_β*U/2)) ) 
Energy[1,1] = 0.0;

#---------------INITIALIZE FOR bicg---------------------------------

	u_r = zeros(Complex{Float64},SpaceSize);
    
 indext = 2;

   for idx = 1:N_t
        
 
        # Projection done
 #	TID = time()
     	 ϕ=Update_Quad(Simplices,Nodes,Mesh2Space,N_of_Simplices,U,ϕ);
#	println("UPDATE ", time()-TID)

 #	TID = time()
        M_ϕ = Reinit_M(Simplices,Nodes,Mesh2Space, M_ϕ, N_of_Simplices);
#	println("Reinit_M ", time()-TID)

# 	TID = time()
        M_ϕ = AssembleMatrix_M_RE(Simplices,Nodes,Mesh2Space,ϕ, M_ϕ, N_of_Simplices,Quad);
#	println("AssembleMatrix_M_RE ", time()-TID)
        
        Matris = M+1im*tau/2*(ϵ*A+β*M_ϕ+M_v)
        RHS =M*U-1im*tau/2*(ϵ*A*U+β*M_ϕ*U+M_v*U);

        println("SOLVING GMRES")
        TID = time()
        r = RHS-Matris*U;
        bicgstabl!(u_r, Matris, r[:],tol=10.0^-14);
        U+=u_r;
        println("TIME ",time()-TID)

        println(" ")

    
#-----------------------------------------------------------
        t = idx*tau
        
        if(T/20*(indext-1)<= t)
  		M_β=Reinit_M(Simplices,Nodes,Mesh2Space, M_β, N_of_Simplices);
		M_β = AssembleMatrix_M_formatFEM(Simplices,Nodes,Mesh2Space,  U[:],M_β, N_of_Simplices,Quad ); 
		Energy[indext,2] = real( dot(U,(0.5*A*U+M_v*U+β*M_β*U/2)) ) 
		Energy[indext,1] = t;

            	indext = indext+1;
       end


    end

Tid = time()-Tid
    
    

 	save(pathway*"RE_beta_"*string(β)*"_Nt_"*string(N_t)*"_Mesh_"*string(Mesh)*".jld","U",U,"Energy",Energy,"Mesh2Space",Mesh2Space,"beta",β,"Tid",Tid)
    
end
