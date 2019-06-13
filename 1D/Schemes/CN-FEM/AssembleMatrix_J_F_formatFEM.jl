function AssembleMatrix_J_F_formatFEM( Simplexes::Array{Int64,2},Nodes::Array{Float64,1},Mesh2Space::Array{Int64,1} ,U_it::Array{Complex{Float64},1},U0::Array{Complex{Float64},1}, 
M_1::SparseMatrixCSC{Float64,Int64},
M_2::SparseMatrixCSC{Float64,Int64},
M_3::SparseMatrixCSC{Float64,Int64},
M_4::SparseMatrixCSC{Float64,Int64}, 
N_of_Simplexes::Int64,
Quad::Dict{String,Array{Float64,1}})
#Calculates <|phi|^2 phi_j, phi_i>
    
    
N_gp=( Int(Quad["N_gp"][1]))


    
 U(X::Float64,C1::Complex{Float64},C2::Complex{Float64}) = (1-X)*C1 +C2*X       

    phi = [0.0 0.0];
    u = 0.0+0.0im;
    u0 = 0.0+0.0im;

#THIS IS STUPID, REPLACE SOMETIME
    N = length(Nodes);
    #translate solution points to whole mesh 
    Data = zeros(N)*0im 
    Data_0 = zeros(N)*0im
    
    loc_1 = zeros(2,2);
    loc_2 = zeros(2,2);
    loc_3 = zeros(2,2);
    loc_4 = zeros(2,2);


    for i = 1:N
       if( Mesh2Space[i] != 0)
          Data[i] = U_it[Mesh2Space[i]]; 
	  Data_0[i] = U0[Mesh2Space[i]];
       end

    end


    for i = 1: N_of_Simplexes

        Simplex= Simplexes[:,i];

        Vertexes = Nodes[Simplex];
        J = Vertexes[2]-Vertexes[1];

        DetJ = abs(J);



        C1 = Data[Simplex[1]];
        C2 = Data[Simplex[2]];
  
        C1_0 = Data_0[Simplex[1]];
        C2_0 = Data_0[Simplex[2]];




        fill!(loc_1,0);
        fill!(loc_2,0);
        fill!(loc_3,0);
        fill!(loc_4,0);

        for gp = 1:N_gp
            X = Quad["X"][gp]; ω =Quad["ω"][gp];
            
	    phi = Φ(X,phi);
	    u = U(X,C1,C2);
            u0 = U(X,C1_0,C2_0);
            for i = 1:2
                for j = 1:2
                    loc_1[i,j] +=  (2*real(u)*(real(u)+real(u0)) + abs(u)^2+abs(u0)^2)*phi[i]*phi[j]*ω;
                
                    loc_2[i,j] += 2*imag(u)*(real(u)+real(u0))*phi[i]*phi[j]*ω;

                    loc_3[i,j] += 2*real(u)*(imag(u)+imag(u0))*phi[i]*phi[j]*ω;

                    loc_4[i,j] += (2*imag(u)*(imag(u)+imag(u0)) + abs(u)^2+abs(u0)^2)*phi[i]*phi[j]*ω;
                end
            end
        end
            loc_1 = loc_1*DetJ;
            loc_2 = loc_2*DetJ;
            loc_3 = loc_3*DetJ;
            loc_4 = loc_4*DetJ;

    
         #Add to Global Matrix
         #Add to Global Matrix
         idx = Mesh2Space[Simplex];


            for i = 1:2
                for j = 1:2
                    if( (idx[i]!=0) & (idx[j]!=0) )  #Om ej kant
                         M_1[idx[i],idx[j]] +=  loc_1[i,j];  
			 M_2[idx[i],idx[j]] +=  loc_2[i,j];       
			 M_3[idx[i],idx[j]] +=  loc_3[i,j];
			 M_4[idx[i],idx[j]] +=  loc_4[i,j];

                     end
                end

            end


    end



        return M_1,M_2,M_3,M_4;

end

