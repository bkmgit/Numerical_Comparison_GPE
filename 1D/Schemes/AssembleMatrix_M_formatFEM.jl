function AssembleMatrix_M_formatFEM( Simplexes::Array{Int64,2},Nodes::Array{Float64,1},Mesh2Space::Array{Int64,1} ,weight::Array{Complex{Float64},1}, M::SparseMatrixCSC{Float64,Int64}, N_of_Simplexes::Int64,Quad::Dict{String,Array{Float64,1}})
#Calculates <|phi|^2 phi_j, phi_i>
    
    N_gp = Int(Quad["N_gp"][1])
#Help functions

 U(X::Float64,C1::Complex{Float64},C2::Complex{Float64}) = (1-X)*C1 +C2*X         
#---------------------------

#initilise variables -------
    phi = [0.0 0.0];
    u = 0.0+0.0im;

    N = length(Nodes);
    #translate solution points to whole mesh 
    Data = zeros(Complex{Float64},N);

    
    loc = zeros(2,2);

#--------------------------------

    for i = 1:N
       if( Mesh2Space[i] != 0)
          Data[i] = weight[Mesh2Space[i]]; 
       end

    end


    for i = 1: N_of_Simplexes

        Simplex = Simplexes[:,i];

        Vertexes = Nodes[Simplex];
         
        J = Vertexes[2]-Vertexes[1];


        DetJ = abs(J);

       

        C1 = Data[Simplex[1]];
        C2 = Data[Simplex[2]];

        fill!(loc,0);

        for gp = 1:N_gp
            X = Quad["X"][gp]; ω =Quad["ω"][gp];

            phi = Φ(X,phi);
            w = abs(U(X,C1,C2))^2;
            for i = 1:2
                for j = 1:2
                    loc[i,j] +=  w*phi[i]*phi[j]*ω;
                end
            end
        end
            loc = loc*DetJ;

    
         #Add to Global Matrix
         #Add to Global Matrix
         idx = Mesh2Space[Simplex];


            for i = 1:2
                for j = 1:2
                    if( (idx[i]!=0) & (idx[j]!=0) )  #Om ej kant
                         M[idx[i],idx[j]] =  M[idx[i],idx[j]]+loc[i,j];       

                     end
                end

            end


    end



        return M   

end

