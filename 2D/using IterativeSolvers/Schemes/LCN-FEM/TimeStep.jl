function TimeStep(Simplices::Array{Int64,2},Nodes::Array{Float64,2},Mesh2Space::Array{Int64,1},N_of_Simplices::Int64,N::Int64,beta::Float64,tau::Float64,A::SparseMatrixCSC{Float64,Int64},M_β::SparseMatrixCSC{Float64,Int64},M::SparseMatrixCSC{Float64,Int64},M_v::SparseMatrixCSC{Float64,Int64},Uhat::Array{Complex{Float64},1},U0::Array{Complex{Float64},1},Energy::Array{Float64,2},ϵ::Float64,Quad::Dict{String,Array{Float64,1}},v::Function,M_L,pathway)

indext = 2;

u_r = zeros(SpaceSize)*0im;
U = Uhat;
    for idx = 1:N_t

    	t = idx*tau;
        
  	M_β=Reinit_M(Simplices,Nodes,Mesh2Space, M_β, N_of_Simplices);
        M_β = AssembleMatrix_M_formatFEM(Simplices,Nodes,Mesh2Space,Uhat,M_β, N_of_Simplices,Quad);


        Matris = -2.0*M+tau*1im*(-ϵ*A-beta*M_β-M_v+Omega*M_L);

        b = -2.0*(M*U0); 
	println("Solving")
	TID = time()
     #   U = Matris\b;

	r = b-Matris*U;
        bicgstabl!(u_r, Matris, r[:],tol=10.0^-14);
        U+=u_r;
	println(time()-TID)

        U = 2.0*U-U0;
        

        Uhat = 0.5*(3.0*U-U0);

        U0 = U;

      if(T/20*(indext-1)<= t)
		M_β=Reinit_M(Simplices,Nodes,Mesh2Space, M_β, N_of_Simplices);
		M_β = AssembleMatrix_M_formatFEM(Simplices,Nodes,Mesh2Space,  U[:],M_β, N_of_Simplices,Quad ); 
		Energy[indext,2] = real( dot(U,(ϵ*A*U+M_v*U+β*M_β*U/2)) ) 
		Energy[indext,1] = t;
		indext +=1;
       end


    


    end
    return U, Energy;
end
