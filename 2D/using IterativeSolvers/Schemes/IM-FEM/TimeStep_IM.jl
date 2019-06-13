function TimeStep_IM(Simplexes::Array{Int64,2},Nodes::Array{Float64,2},Mesh2Space::Array{Int64,1},N_of_Simplexes::Int64,
        N::Int64,
        β::Float64,
	N_t,
        tau::Float64,
        A::SparseMatrixCSC{Float64,Int64},
        M::SparseMatrixCSC{Float64,Int64},
        M_v::SparseMatrixCSC{Float64,Int64},
        U0::Array{Complex{Float64},1},
        Energy::Array{Float64,2},
        ϵ::Float64,
        T::Float64,
        θ::Float64,
        SpaceSize::Int64,
        Quad::Dict{String,Array{Float64,1}},
        Mesh::Int64)
    
max_it =8;
    
diag1 = 0; diag2 = 1; diag3 = Mesh-2; diag4 = Mesh-1;
vec1 = ones(SpaceSize)*eps(1/1000000);  vec2=vec1[1:SpaceSize-diag2]; vec3 = vec1[1:SpaceSize-diag3]; vec4 = vec1[1:SpaceSize-diag4];

∂F_R∂ξ_R=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;

∂F_R∂ξ_I=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;

∂F_I∂ξ_R=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;

∂F_I∂ξ_I=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;

M_β = spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;


F_real = zeros(SpaceSize);
fill!(F_real,0);
F_imag = zeros(SpaceSize)
fill!(F_imag,0);
    
U = U0;

#----------------------------Measure Initial Energy ---------------------------------------
M_β = AssembleMatrix_M_formatFEM(Simplexes,Nodes,Mesh2Space,  U[:],M_β, N_of_Simplexes,Quad );
Energy[1,2] = real((U'*(A+M_v+β*M_β /2)*U));
indext = 2;
#------------------------------------------------------------------------------------------

    for idx = 1:N_t
	
	t = tau*idx;

          g  = [-M*imag(U0)+tau*(1-θ)*(A*real(U0)+M_v*real(U0));
                 M*real(U0)+tau*(1-θ)*(A*imag(U0)+M_v*imag(U0))]; # CONSTANT DURING NEWTON STEP, NON linearity unfortunately included in F

        #Newton Step
        ΔU = ones(2*SpaceSize);
	
        it=1;

        while( ((norm(ΔU)/norm(U))>10.0^-8) & (it<max_it))
            it+=1;
            ∂F_R∂ξ_R, ∂F_R∂ξ_I, ∂F_I∂ξ_R, ∂F_I∂ξ_I = AssembleMatrix_J_F_formatFEM_IM(Simplexes,Nodes,Mesh2Space,U,U0, ∂F_R∂ξ_R,∂F_R∂ξ_I,∂F_I∂ξ_R,∂F_I∂ξ_I, N_of_Simplexes,θ,Quad); 

            F_real,F_imag = Assemble_F_formatFEM_IM(Simplexes,Nodes,Mesh2Space ,U,U0, F_real,F_imag, N_of_Simplexes,θ,Quad)


            G = [-M*imag(U)-tau*(θ*A*real(U)+θ*M_v*real(U)+β*F_real);
                  M*real(U)-tau*(θ*A*imag(U)+θ*M_v*imag(U)+β*F_imag)];
            
            J_G = [-tau*(θ*A+θ*M_v+β*∂F_R∂ξ_R) -M-tau*β*∂F_R∂ξ_I ; 
                    M-tau*β*∂F_I∂ξ_R           -tau*(θ*A+θ*M_v+β*∂F_I∂ξ_I)];


            sol = [real(U);imag(U)];
 
      #     println("solving GMRES")
        #    tic()
            r = G-g;
            bicgstabl!(ΔU, J_G, r[:],tol=10.0^-15);
            
        #    toc()
   
            val = 0.0;
	   
            sol = sol - ΔU;

            U = sol[1:SpaceSize]+1im*sol[SpaceSize+1:end];

            fill!(F_real,0); fill!(F_imag,0); 
            

           ∂F_R∂ξ_R= Reinit_M(Simplexes,Nodes,Mesh2Space, ∂F_R∂ξ_R, N_of_Simplexes);
           ∂F_R∂ξ_I= Reinit_M(Simplexes,Nodes,Mesh2Space, ∂F_R∂ξ_I, N_of_Simplexes);
           ∂F_I∂ξ_R= Reinit_M(Simplexes,Nodes,Mesh2Space, ∂F_I∂ξ_R, N_of_Simplexes);
           ∂F_I∂ξ_I= Reinit_M(Simplexes,Nodes,Mesh2Space, ∂F_I∂ξ_I, N_of_Simplexes);
           

        end
         

        if(T/20*(indext-1)<= t)
  		M_β=Reinit_M(Simplices,Nodes,Mesh2Space, M_β, N_of_Simplices);
		M_β = AssembleMatrix_M_formatFEM(Simplices,Nodes,Mesh2Space,  U[:],M_β, N_of_Simplices,Quad ); 
		Energy[indext,2] = real( dot(U,(A*U+M_v*U+β*M_β*U/2)) ) 
		Energy[indext,1] = t;

            	indext = indext+1;
       end
 
            U0=U;
println(t)

    end
     
    return U0, Energy;
end


