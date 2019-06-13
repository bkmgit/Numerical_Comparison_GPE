function Relax(Simplexes,Nodes,Mesh2Space,SpaceSize,N_of_Simplexes,N_of_Nodes,T,U0,tau,Mesh,pathway::String,v::Function)

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

ϕ = zeros(6,N_of_Simplexes);
#Phiold = zeros(6,N_of_Simplexes);
ϕ=Start_Quad(Simplexes,Nodes,Mesh2Space,N_of_Simplexes,U0,ϕ);
σ = 1;
    
Energy = zeros(21,2);
Energy[1,1] = E_Norm( Simplexes,Nodes,Mesh2Space,U0[:] ,N_of_Nodes,β,N_of_Simplexes,ϵ,v,Quad);
Energy[1,2] = 0.0;
    
#-----------------------------Assemble Matrices ------------------------------------------------------------------- 

    func(x) =1.0;

#tic()
    M = AssembleMatrix_M(Simplexes,Nodes,Mesh2Space, func,M,N_of_Simplexes,Quad);
#toc()

#tic()
    A=AssembleMatrix_A(Simplexes,Nodes,Mesh2Space,A, N_of_Simplexes);
#toc()

#tic()
    M_v = AssembleMatrix_M(Simplexes,Nodes,Mesh2Space,v,M_v, N_of_Simplexes,Quad);
#toc()
    M_β = AssembleMatrix_M_formatFEM(Simplexes,Nodes,Mesh2Space,  U0[:],M_β, N_of_Simplexes,Quad ); 
    w =abs.(U0).^(2*σ);
    # -----------------------------------------------------------------
    I = spdiagm(0=>ones(SpaceSize));
    p = zeros(SpaceSize,1);
    
    U = U0;
  # Uold = U0;  
    indext=2;
    indext2=2;













    
    #INITIALIZE FOR bicg

u_r = zeros(Complex{Float64},SpaceSize);










#----------------------------TAKE HALF STEP-----------------------------------#


        M_ϕ = Reinit_M(Simplexes,Nodes,Mesh2Space, M_ϕ, N_of_Simplexes);
#	println("Reinit_M ", time()-TID)

# 	TID = time()
        M_ϕ = Assemble_M(Simplexes,Nodes,Mesh2Space,ϕ, M_ϕ, N_of_Simplexes,Quad);
#	println("Assemble_M ", time()-TID)
        
        Matris = M+1im*tau/2/2*(ϵ*A+β*M_ϕ+M_v)
        RHS =M*U-1im*tau/2/2*(ϵ*A*U+β*M_ϕ*U+M_v*U);

	r = RHS-Matris*U0;
        bicgstabl!(u_r, Matris, r[:],tol=10.0^-14);
        U0+=u_r;	

ϕ=Start_Quad(Simplexes,Nodes,Mesh2Space,N_of_Simplexes,U0,ϕ);
#---------------------------------------------------------------------------#    
 
   for t = tau:tau:T
        

        # Projection done
 #	TID = time()
     	 ϕ=Update_Quad(Simplexes,Nodes,Mesh2Space,N_of_Simplexes,U,ϕ);
#	println("UPDATE ", time()-TID)

 #	TID = time()
        M_ϕ = Reinit_M(Simplexes,Nodes,Mesh2Space, M_ϕ, N_of_Simplexes);
#	println("Reinit_M ", time()-TID)

# 	TID = time()
        M_ϕ = Assemble_M(Simplexes,Nodes,Mesh2Space,ϕ, M_ϕ, N_of_Simplexes,Quad);
#	println("Assemble_M ", time()-TID)
        
        Matris = M+1im*tau/2*(ϵ*A+β*M_ϕ+M_v)
        RHS =M*U-1im*tau/2*(ϵ*A*U+β*M_ϕ*U+M_v*U);

        println("SOLVING GMRES")
        TID = time()
        r = RHS-Matris*U;
        bicgstabl!(u_r, Matris, r[:],tol=10.0^-14);
        U+=u_r;
        println("TIME ",time()-TID)
       # tic()
        # U = Matris\RHS;
       # println(norm(U-U1))
      #  U1 = qrfact(Matris)\RHS;
       # U2 = lufact(Matris)\RHS;
        #toc()
        println(" ")

       #Mass = (U'*(M*U));
	#println(real(Mass[1]))
        
#println(norm(U1-U))
     

       #-----------------------------------------------------------
        
        
        if(T/20*(indext-1)<= t)
#		        TID = time()
  		     M_β=Reinit_M(Simplexes,Nodes,Mesh2Space, M_β, N_of_Simplexes);
		     M_β = AssembleMatrix_M_formatFEM(Simplexes,Nodes,Mesh2Space,  U[:],M_β, N_of_Simplexes,Quad ); 
		     Energy[indext,1] = real((U'*(0.5*A+M_v+2300*M_β /2)*U))
            #Energy[indext,1] = E_Norm( Simplexes,Nodes,Mesh2Space,U[:] ,N_of_Nodes,β,N_of_Simplexes,ϵ,v,Quad);
#		        println("Energy ", time()-TID)
            Energy[indext,2] = t;
         #   E_Norm_Relax(Simplexes,Nodes,Mesh2Space, Uold[:] ,N_of_Nodes,N_of_Simplexes , ϵ,v,Phiold,ϕ,Quad)

            
          indext = indext+1;
       end


	#if(T/1000*(indext2-1)<= t)
		#println(maximum(abs.(U).^2))
           # save(pathway*"Sol"*string(indext2)* ".jld", "U",U)
       #     npzwrite("./Test/Sol"*string(indext2)*".npy", abs.(U).^2)
       #     indext2 = indext2+1;
       #	end
    


       # Phiold[:,:] = ϕ[:,:];
       # Uold[:]= U[:];
    end

Tid -=time();
    
    

 	save(pathway*"Relaxation_beta"*string(β)*"_tau"*string(tau)*"_Mesh"*string(Mesh)*".jld","U",U,"Energy",Energy,"Mesh2Space",Mesh2Space,"beta",β,"Tid",Tid)
   # return U;
    
    
    
    
    
end
