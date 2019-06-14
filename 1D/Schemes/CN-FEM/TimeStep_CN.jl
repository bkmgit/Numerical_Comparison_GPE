function TimeStep_CN(Simplices::Array{Int64,2},Nodes::Array{Float64,1},Mesh2Space::Array{Int64,1},N_of_Simplices::Int64,
        N::Int64,
        tau::Float64,
	N_t,
        A::SparseMatrixCSC{Float64,Int64},
        M::SparseMatrixCSC{Float64,Int64},
        M_v::SparseMatrixCSC{Float64,Int64},
        U0::Array{Complex{Float64},1},
        Energy::Array{Float64,2},
        SpaceSize::Int64,
        Quad::Dict{String,Array{Float64,1}},
	v::Function,
	E::Array{Float64,1})

indext = 2;
max_it = 8; #Maximum number of newton iterations

#---------------------Initiliaze matrices for Newton step---------------------
diag1 = 0; diag2 = 1; 
vec1 = ones(SpaceSize)*eps(1/10000);  vec2=vec1[1:SpaceSize-diag2]; 

∂F_R∂ξ_R = spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);
;
∂F_R∂ξ_I = spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);
;
∂F_I∂ξ_R = spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);
;
∂F_I∂ξ_I = spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);

M_b = spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);

F_real = zeros(SpaceSize);
fill!(F_real,0);
F_imag = zeros(SpaceSize)
fill!(F_imag,0);

#----------------------------------------------------------------------------

U = U0;
idx_E = 1;

t = 0.

    for idx_t = 1:N_t



	#---------------------Energy Measurement ---------------------
       if( idx_t == E[idx_E]) 
	    M_b=Reinit_M(Simplices,Nodes,Mesh2Space, M_b, N_of_Simplices);
	    M_b = AssembleMatrix_M_formatFEM(Simplices,Nodes,Mesh2Space,U, M_b, N_of_Simplices,Quad)
	    Energy[idx_E,2] = real(U'*(A*U)+β*U'*(M_b*U)/2.0+U'*M_v*U)
 	    Energy[idx_E,1] = t
	    idx_E +=1;	 
       end
	#---------------------------------------------------------------
        

    g  = [-M*imag(U0)+tau/2*(A*real(U0)+M_v*real(U0) ); M*real(U0) + tau/2*(A*imag(U0)+M_v*imag(U0))]; # CONSTANT DURING NEWTON STEP

        #Newton Step

	# SOME NOTES
	#i<u,φ>-τ/2(<∇u,∇φ> + <Vu,φ> + <(|u|^2+|u0|^2)*(u+u0),φ> = i<u0,φ>+ +τ/2(<∇u0,∇φ>+<Vu0,φ>) is solved as:
	#G(u) = g(u0) where G(u) = LHS = i<u,φ>-τ/2(<∇u,∇φ> + <Vu,φ> + <(|u|^2+|u0|^2)*(u+u0),φ>
	#and g(u0) = RHS
	# u = u- J_G\(G(u)-g(u0))

        ΔU = 1;
        it=0;
        while( (norm(ΔU)>10.0^-10) & (it<max_it))
            it+=1;
            ∂F_R∂ξ_R, ∂F_R∂ξ_I, ∂F_I∂ξ_R, ∂F_I∂ξ_I = AssembleMatrix_J_F_formatFEM(Simplices,Nodes,Mesh2Space,U,U0, ∂F_R∂ξ_R,∂F_R∂ξ_I,∂F_I∂ξ_R,∂F_I∂ξ_I, N_of_Simplices,Quad); 

            F_real,F_imag = Assemble_F_formatFEM(Simplices,Nodes,Mesh2Space ,U,U0, F_real,F_imag, N_of_Simplices,Quad)

            G = [-M*imag(U)-tau/2*(A*real(U)+M_v*real(U)+β/2*F_real); M*real(U)-tau/2*(A*imag(U)+M_v*imag(U) + β/2*F_imag)];
            J_G = [-tau/2*(A+M_v+β/2*∂F_R∂ξ_R) -M-tau/2*β/2*∂F_R∂ξ_I ; M-tau/2*β/2*∂F_I∂ξ_R -tau/2*(A+M_v+β/2*∂F_I∂ξ_I)];


            sol = [real(U);imag(U)];
            ΔU = J_G\(G-g);
  
            sol = sol - ΔU;

            U = sol[1:SpaceSize]+1im*sol[SpaceSize+1:end];
	    t = tau*idx_t;

           F_real *= 0; F_imag *= 0;
 
           ∂F_R∂ξ_R= Reinit_M(Simplices,Nodes,Mesh2Space, ∂F_R∂ξ_R, N_of_Simplices);
           ∂F_R∂ξ_I= Reinit_M(Simplices,Nodes,Mesh2Space, ∂F_R∂ξ_I, N_of_Simplices);
           ∂F_I∂ξ_R= Reinit_M(Simplices,Nodes,Mesh2Space, ∂F_I∂ξ_R, N_of_Simplices);
           ∂F_I∂ξ_I= Reinit_M(Simplices,Nodes,Mesh2Space, ∂F_I∂ξ_I, N_of_Simplices);
           

        end


        t = idx_t*tau;

      	U0=U;
   end


	#Measure Final Energy
	M_b=Reinit_M(Simplices,Nodes,Mesh2Space, M_b, N_of_Simplices);
	M_b = AssembleMatrix_M_formatFEM(Simplices,Nodes,Mesh2Space,U, M_b, N_of_Simplices,Quad)
	Energy[idx_E,2] = real(ϵ*U'*(A*U)+β*U'*(M_b*U)/2.0+U'*M_v*U)
	Energy[idx_E,1] = T
     
    return U0, Energy;
end


