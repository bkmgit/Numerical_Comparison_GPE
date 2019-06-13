function TimeStep_CN(Simplices::Array{Int64,2},Nodes::Array{Float64,2},Mesh2Space::Array{Int64,1},N_of_Simplices::Int64,
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
        SpaceSize::Int64,
        Quad::Dict{String,Array{Float64,1}},
        Mesh::Int64)

indext = 2;
    
Max_it = 8;


diag1 = 0; diag2 = 1; diag3 = Mesh-2; diag4 = Mesh-1;
vec1 = ones(SpaceSize)*eps(1/1000000);  vec2=vec1[1:SpaceSize-diag2]; vec3 = vec1[1:SpaceSize-diag3]; vec4 = vec1[1:SpaceSize-diag4];

∂F_R∂ξ_R=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;

∂F_R∂ξ_I=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;

∂F_I∂ξ_R=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;

∂F_I∂ξ_I=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;

M_β=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;




F_real = zeros(SpaceSize);
F_imag = zeros(SpaceSize)

    
U = U0;

#----------------------------Measure Initial Energy ---------------------------------------
M_β = AssembleMatrix_M_formatFEM(Simplexes,Nodes,Mesh2Space,  U[:],M_β, N_of_Simplexes,Quad );
Energy[1,2] = real((U'*(0.5*A+M_v+β*M_β /2)*U));
indext = 2;
#------------------------------------------------------------------------------------------


    for idx = 1:N_t
	t = tau*idx
        

    g  = [-M*imag(U0)+tau/2*(A*real(U0)+M_v*real(U0) ); 
           M*real(U0)+tau/2*(A*imag(U0)+M_v*imag(U0))]; # CONSTANT DURING NEWTON STEP, NON linearity unfortunately included in F

        #Newton Step
        ΔU = ones(2*SpaceSize);
	it = 1;
        while(norm(ΔU)>10.0^-8 && it <8)

           it+=1;

            ∂F_R∂ξ_R, ∂F_R∂ξ_I, ∂F_I∂ξ_R, ∂F_I∂ξ_I = AssembleMatrix_J_F_formatFEM(Simplices,Nodes,Mesh2Space,U,U0, ∂F_R∂ξ_R,∂F_R∂ξ_I,∂F_I∂ξ_R,∂F_I∂ξ_I, N_of_Simplices,Quad); 

            F_real,F_imag = Assemble_F_formatFEM(Simplices,Nodes,Mesh2Space ,U,U0, F_real,F_imag, N_of_Simplices,Quad)

            G = [-M*imag(U)-tau/2*(A*real(U)+M_v*real(U)+β/2*F_real);
                  M*real(U)-tau/2*(A*imag(U)+M_v*imag(U)+β/2*F_imag)];
            
            J_G = [-tau/2*(A+M_v+β/2*∂F_R∂ξ_R) -M-tau/2*β/2*∂F_R∂ξ_I ; 
                   M-tau/2*β/2*∂F_I∂ξ_R        -tau/2*(A+M_v+β/2*∂F_I∂ξ_I)];


            sol = [real(U);imag(U)];

#            ΔU = J_G\(G-g);

     #	    println("solving GMRES")
            r = G-g;
            bicgstabl!(ΔU, J_G, r[:],tol=10.0^-15);

            sol = sol - ΔU;

            U = sol[1:SpaceSize]+1im*sol[SpaceSize+1:end];

            F_real *= 0; F_imag *= 0; #∂F_R∂ξ_R *= eps(1/100000);∂F_R∂ξ_I *= eps(1/100000); ∂F_I∂ξ_R *= eps(1/100000) ; ∂F_I∂ξ_I *= eps(1/100000); 
	   ∂F_R∂ξ_R= Reinit_M(Simplices,Nodes,Mesh2Space, ∂F_R∂ξ_R, N_of_Simplices);
           ∂F_R∂ξ_I= Reinit_M(Simplices,Nodes,Mesh2Space, ∂F_R∂ξ_I, N_of_Simplices);
           ∂F_I∂ξ_R= Reinit_M(Simplices,Nodes,Mesh2Space, ∂F_I∂ξ_R, N_of_Simplices);
           ∂F_I∂ξ_I= Reinit_M(Simplices,Nodes,Mesh2Space, ∂F_I∂ξ_I, N_of_Simplices);
           

        end
        	println(t)

        
        if(T/20*(indext-1)<= t)
  		M_β=Reinit_M(Simplices,Nodes,Mesh2Space, M_β, N_of_Simplices);
		M_β = AssembleMatrix_M_formatFEM(Simplices,Nodes,Mesh2Space,  U[:],M_β, N_of_Simplices,Quad ); 
		Energy[indext,2] = real( dot(U,(A*U+M_v*U+β*M_β*U/2)) ) 
		Energy[indext,1] = t;

            	indext = indext+1;
       end
            U0=U;
    end
     
    return U0, Energy;
end






# SOME NOTES
#i<u,φ>-τ/2(<∇u,∇φ> + <Vu,φ> + <(|u|^2+|u0|^2)*(u+u0),φ> = i<u0,φ>+ +τ/2(<∇u0,∇φ>+<Vu0,φ>) is solved as:
#G(u) = g(u0) where G(u) = LHS = i<u,φ>-τ/2(<∇u,∇φ> + <Vu,φ> + <(|u|^2+|u0|^2)*(u+u0),φ>
#and g(u0) = RHS
# u = u- J_G\(G(u)-g(u0))
