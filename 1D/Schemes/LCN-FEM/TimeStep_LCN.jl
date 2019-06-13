function TimeStep_LCN(Simplices::Array{Int64,2},
Nodes::Array{Float64,1},
Mesh2Space::Array{Int64,1},
N_of_Simplices::Int64,
N::Int64,
beta::Float64,
tau::Float64,
N_t,
A::SparseMatrixCSC{Float64,Int64},
M_β::SparseMatrixCSC{Float64,Int64},
M::SparseMatrixCSC{Float64,Int64},
M_v::SparseMatrixCSC{Float64,Int64},
Uhat::Array{Complex{Float64},1},
U0::Array{Complex{Float64},1},
ϵ::Float64,
Quad::Dict{String,Array{Float64,1}},
v::Function)


U = U0;
t = 0.;

diag1 = 0; diag2 = 1;
vec1 = ones(SpaceSize)*eps(1/10000);  vec2=vec1[1:SpaceSize-diag2];
M_βE = spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);

#-----------------------------Do 200 Energy Measurements------------------------------------------"
N_E = 200;
if(N_t>N_E)
	E =[1+N_t/N_E*i for i =0:N_E]
	E = round.(E)
	Energy = zeros(N_E+1,2);
else
	E = [1+i for i = 0:N_t]
	E = Array{Int}(round.(E))
	Energy = zeros(Int(N_t)+1,2);
end 

idx_E = 1
#-----------------------------------------------------------------------------------------------------

    for idx_t = 1:N_t
	
	


	M_β=Reinit_M(Simplices,Nodes,Mesh2Space, M_β, N_of_Simplices);
        M_β = AssembleMatrix_M_formatFEM(Simplices,Nodes,Mesh2Space,Uhat[:],M_β, N_of_Simplices,Quad);


	# Measure Energy
	if( idx_t == E[idx_E]) 
	    	M_βE=Reinit_M(Simplices,Nodes,Mesh2Space, M_βE, N_of_Simplices);
	        M_βE = AssembleMatrix_M_formatFEM(Simplices,Nodes,Mesh2Space,U[:],M_βE, N_of_Simplices,Quad);
		Energy[idx_E,2] = real(ϵ*U'*(A*U)+β*U'*(M_β*U)/2.0+U'*M_v*U)
		Energy[idx_E,1] = t
	        idx_E +=1;
	end


        Matris = -2.0*M+1im*tau*(-ϵ*A-beta*M_β-M_v);

        b = -2.0*(M*U0); 


        U = Tridiagonal(Matris)\b;


        U = 2.0*U-U0;
        
	
        Uhat = 0.5*(3.0*U-U0);

        U0 = U;
	
	t = idx_t*tau;

   


    end


	#Measure Final Energy
	Energy[idx_E,2] = real(ϵ*U'*(A*U)+β*U'*(M_β*U)/2.0+U'*M_v*U)
	Energy[idx_E,1] = t


    return U0, Energy;
end
